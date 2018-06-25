close all
clear
clc

% load('C:\Data\2018_05_keyboard\Tests\20180604T161434_longTraining.mat')
load('fileList.mat');
obj=MI_MIkeyboard.joinSessions(f);

% Recover correct selCounter from UDP transmissions
selCounter=cellfun(@(x)str2double(char(x)),obj.outputLog.keyUDPlogOut);

% Subsample selCounter and set failed trials counter to 0
flexPoints=find((selCounter(1:end-2)-selCounter(2:end-1)).*(selCounter(3:end)-selCounter(2:end-1))>=0)+1;
trialStart=flexPoints(selCounter(flexPoints)<selCounter(flexPoints+1));
trialEnd=flexPoints(selCounter(flexPoints)>selCounter(flexPoints-1));
trialLbls=(selCounter(trialEnd)==1)';
negTrials=find(trialLbls==0);
for currTrial=1:length(negTrials)
    selCounter(trialStart(negTrials(currTrial)):trialEnd(negTrials(currTrial)))=0;
end
selCounterShort=selCounter(1:2:end);
selCounterShortTimes=obj.outputLog.time(1:2:end);

% Divide collected data in different trials
artifactLimit=15;
timeWindowLength=.2;
selCounterShort(selCounterShortTimes<=artifactLimit)=[];
selCounterShortTimes(selCounterShortTimes<=artifactLimit)=[];
trialData=cell(length(selCounterShort),1);
for currTrial=1:length(selCounterShort)
    relLims=[selCounterShortTimes(currTrial)-timeWindowLength,selCounterShortTimes(currTrial)];
    relIdx=(obj.rawData.time>relLims(1))&(obj.rawData.time<=relLims(2));
    relData=obj.rawData.data(relIdx,:)';
    trialData{currTrial}=relData;
end

% Create network
nNeurons= 40;
maxEpochs = 100;
validationChecks=15;
miniBatchSize = round(length(selCounterShort)/5);
inputSize = size(trialData{1},1);
trainErr=zeros(length(nNeurons),1); %#ok<*PREALL>
testErr=zeros(length(nNeurons),1);

% Define network structure and training options
layers=[ ...
    sequenceInputLayer(inputSize)
    lstmLayer(nNeurons,'OutputMode','last')
    fullyConnectedLayer(10)
    dropoutLayer(0.5)
    fullyConnectedLayer(1)
    regressionLayer];
options = trainingOptions('adam', ...
    'MaxEpochs',1, ...
    'MiniBatchSize',miniBatchSize, ...
    'InitialLearnRate',0.01, ...
    'GradientThreshold',1, ...
    'Shuffle','never', ...
    'Verbose',0, ...
    'ExecutionEnvironment','gpu');

% Define blocks for k-fold validation
C.NumTestSets=5;
% Random grouping
% C.groups=ceil(rand(size(selCounterShort))*C.NumTestSets);
% Subsequent grouping
C.groups=ceil(linspace(1/length(selCounterShort),C.NumTestSets,length(selCounterShort)));
C.training=@(currGroup)C.groups~=currGroup;
C.test=@(currGroup)C.groups==currGroup;
selCounterShortEst=zeros(size(selCounterShort));

computeErr=@(x,y)mean((x-y).^2);
for currFold=1:C.NumTestSets
    % Recover training and testing sets
    trainData=trialData(C.training(currFold));
    trainLbls=selCounterShort(C.training(currFold));
    testData=trialData(C.test(currFold));
    testLbls=selCounterShort(C.test(currFold));
    
    % Split training set in training and validation
    valSetProportion=.5;
    rndOrdr=randperm(length(trainData));
    valIdx=rndOrdr(1:floor(length(trainData)*valSetProportion));
    trainIdx=setdiff(1:length(trainLbls),valIdx);
    
    valData=trainData(valIdx);
    valLbls=trainLbls(valIdx);
    trainData=trainData(trainIdx);
    trainLbls=trainLbls(trainIdx);
    
    % Train network
    currEpoch=1;
    valErr=zeros(maxEpochs,1);
    fltrdValErr=ones(maxEpochs,1)*Inf;
    while true
        if exist('net','var')
            net=trainNetwork(trainData',trainLbls',net.Layers,options);
        else
            net=trainNetwork(trainData',trainLbls',layers,options);
        end
        
        % Generate validation lbls
        valLblsEst=predict(net,valData)';
        trainLblsEst=predict(net,trainData)';
        valErr(currEpoch)=computeErr(valLbls,valLblsEst);
        fltrdValErr(currEpoch)=median(valErr(max(1,currEpoch-4):currEpoch));
        if fltrdValErr(currEpoch)==min(fltrdValErr)
            bestNet=net;
        end
        trainErr=computeErr(trainLbls,trainLblsEst);
        fprintf('Epoch: %d/%d, train Err: %0.2f, validation Err: %0.2f\n',currEpoch,maxEpochs,trainErr,min(fltrdValErr));
        
        % Test for ending conditions
        currEpoch=currEpoch+1;
        if ((currEpoch==maxEpochs)||(find(fltrdValErr==min(fltrdValErr),1,'last')+validationChecks<currEpoch))
            break
        end
    end
    net=bestNet;
    
    % Generate test lbls
    trainLblsEst=predict(net,trainData)';
    testLblsEst=predict(net,testData)';
    selCounterShortEst(C.test(currFold))=testLblsEst;
    
    % Compute resulting accuracies on training and testing sets
    trainErr=computeErr(trainLbls,trainLblsEst);
    testErr=computeErr(testLbls,testLblsEst);
    fprintf('Train Err: %0.2f; validation Err: %0.2f; test Err: %0.2f\n\n',trainErr,max(fltrdValErr),testErr);
    clear net
end
