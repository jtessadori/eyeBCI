close all
clear
clc

% load('C:\Data\2018_05_keyboard\Tests\20180604T161434_longTraining.mat')
load('fileList.mat');
obj=MI_MIkeyboard.joinSessions(f);

% Recover correct selCounter from UDP transmissions
selCounter=cellfun(@(x)str2double(char(x)),obj.outputLog.UDPlogOut);

% Recover beginning and end of key selections from selCounter, as well as
% labels
flexPoints=find((selCounter(1:end-2)-selCounter(2:end-1)).*(selCounter(3:end)-selCounter(2:end-1))>0)+1;
flexPoints=union(flexPoints,find(diff(selCounter==0)==-1));
trialStart=obj.outputLog.time(flexPoints(selCounter(flexPoints)<=selCounter(flexPoints+1)));
trialEnd=flexPoints(selCounter(flexPoints)>=selCounter(flexPoints+1));
trialLbls=(selCounter(trialEnd)==1)';
trialEnd=obj.outputLog.time(trialEnd);

% Divide collected data in different trials
artifactLimit=15;
timeWindowLength=.4;
trialLbls(trialStart<=artifactLimit)=[];
trialEnd(trialStart<=artifactLimit)=[];
trialStart(trialStart<=artifactLimit)=[];
trialData=cell(length(trialLbls),1);
for currTrial=1:length(trialLbls)
    relLims=[trialStart(currTrial)-timeWindowLength,trialStart(currTrial)];
%     relLims=[trialEnd(currTrial)-timeWindowLength,trialEnd(currTrial)];
    relIdx=(obj.rawData.time>relLims(1))&(obj.rawData.time<=relLims(2));
    relData=obj.rawData.data(relIdx,:)';
    trialData{currTrial}=relData;
%     for currCh=1:size(relData,1)
%         trialData{currTrial}(currCh,:)=resample(relData(currCh,:),64,obj.fs);
%     end
end

% Create network
nNeurons= 100;
maxEpochs = 100;
validationChecks=15;
miniBatchSize = round(length(trialLbls)/5);
inputSize = size(trialData{1},1);
numClasses = 2;
trainBACC=zeros(length(nNeurons),1); %#ok<*PREALL>
testBACC=zeros(length(nNeurons),1);

% Define network structure and training options
layers=[ ...
    sequenceInputLayer(inputSize)
    lstmLayer(nNeurons,'OutputMode','last')
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];
options = trainingOptions('adam', ...
    'ExecutionEnvironment','auto', ...
    'MaxEpochs',1, ...
    'MiniBatchSize',miniBatchSize, ...
    'GradientThreshold',1, ...
    'Verbose',0);%,...
%'Plots','training-progress');

% Define blocks for k-fold validation
C.NumTestSets=20;
% Random grouping
C.groups=ceil(rand(size(trialLbls))*C.NumTestSets);
% Subsequent grouping
% C.groups=ceil(linspace(1/length(trialLbls),C.NumTestSets,length(trialLbls)));
C.training=@(currGroup)C.groups~=currGroup;
C.test=@(currGroup)C.groups==currGroup;
trialLblsEst=zeros(size(trialLbls));

testAcc=@(x,y)(sum((x==1).*(y==1))./sum(x==1)+sum((x==0).*(y==0))./sum(x==0))*.5;
for currFold=1:C.NumTestSets
    % Recover training and testing sets
    trainData=trialData(C.training(currFold));
    trainLbls=trialLbls(C.training(currFold));
    testData=trialData(C.test(currFold));
    testLbls=trialLbls(C.test(currFold));
    
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
    valBACC=zeros(maxEpochs,1);
    fltrdValBACC=zeros(maxEpochs,1);
    while true
        if exist('net','var')
            net=trainNetwork(trainData',categorical(trainLbls),net.Layers,options);
        else
            net=trainNetwork(trainData',categorical(trainLbls),layers,options);
        end
        
        % Generate validation lbls
        valLblsEst=classify(net,valData);
        trainLblsEst=classify(net,trainData);
        valBACC(currEpoch)=testAcc(valLbls,double(valLblsEst)-1);
        fltrdValBACC(currEpoch)=median(valBACC(max(1,currEpoch-4):currEpoch));
        if fltrdValBACC(currEpoch)==max(fltrdValBACC)
            bestNet=net;
        end
        trainBACC=testAcc(trainLbls,double(trainLblsEst)-1);
        fprintf('Epoch: %d/%d, train BACC: %0.2f, validation BACC: %0.2f\n',currEpoch,maxEpochs,trainBACC,fltrdValBACC(currEpoch));
        
        % Test for ending conditions
        currEpoch=currEpoch+1;
        if ((currEpoch==maxEpochs)||(find(fltrdValBACC==max(fltrdValBACC),1,'last')+validationChecks<currEpoch))
            break
        end
    end
    net=bestNet;
    
    % Generate test lbls
    trainLblsEst=classify(net,trainData);
    testLblsEst=classify(net,testData);
    trialLblsEst(C.test(currFold))=double(testLblsEst)-1;
    
    % Compute resulting accuracies on training and testing sets
    trainBACC=testAcc(trainLbls,double(trainLblsEst)-1);
    testBACC=testAcc(testLbls,double(testLblsEst)-1);
    fprintf('Train BACC: %0.2f; validation BACC: %0.2f; test BACC: %0.2f\n\n',trainBACC,max(fltrdValBACC),testBACC);
    clear net
end
