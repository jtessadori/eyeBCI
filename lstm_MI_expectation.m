close all
clear
clc

% load('C:\Data\2018_05_keyboard\Tests\20180604T161434_longTraining.mat')
load('fileList.mat');
obj=MI_MIkeyboard.joinSessions(f);
% load('C:\Data\2018_06_MIkeyboard\S_Francesco\20180612T162658.mat')
% load('C:\Data\2018_06_MIkeyboard\S_Chang\20180614T102919_training.mat')

% Recover correct selCounter from UDP transmissions
selCounter=cellfun(@(x)str2double(char(x)),obj.outputLog.keyUDPlogOut);

% Recover beginning and end of key selections from selCounter, as well as
% labels
flexPoints=find((selCounter(1:end-2)-selCounter(2:end-1)).*(selCounter(3:end)-selCounter(2:end-1))>=0)+1;
trialStartIdx=flexPoints(selCounter(flexPoints)<selCounter(flexPoints+1));
trialEndIdx=flexPoints(selCounter(flexPoints)>selCounter(flexPoints-1));
trialStartTime=obj.outputLog.time(trialStartIdx);
trialEndTime=obj.outputLog.time(trialEndIdx);
trialLbls=(selCounter(trialEndIdx)==1)';

% Divide collected data in different trials
artifactLimit=15;
timeWindowLength=2;
trialLbls(trialStartTime<=artifactLimit)=[];
trialEndTime(trialStartTime<=artifactLimit)=[];
trialStartTime(trialStartTime<=artifactLimit)=[];
trialData=cell(length(trialLbls),1);
for currTrial=1:length(trialLbls)
    relLims=[trialStartTime(currTrial)-timeWindowLength,trialStartTime(currTrial)];
    relIdx=(obj.rawData.time>relLims(1))&(obj.rawData.time<=relLims(2));
    relData=obj.rawData.data(relIdx,:)';
    trialData{currTrial}=relData;
end
toBeRemoved=cellfun(@(x)isempty(x),trialData);
trialData(toBeRemoved)=[];
trialLbls(toBeRemoved)=[];
trialLbls=double(trialLbls);

% Create network
nNeurons= 40;
maxEpochs = 100;
validationChecks=15;
miniBatchSize = round(length(trialLbls)/5);
inputSize = size(trialData{1},1);
numClasses = 2;
trainBAcc=zeros(length(nNeurons),1); %#ok<*PREALL>
testBAcc=zeros(length(nNeurons),1);

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
C.NumTestSets=5;
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
    valBAcc=zeros(maxEpochs,1);
    fltrdValBAcc=zeros(maxEpochs,1);
    while true
        if exist('net','var')
            net=trainNetwork(trainData',categorical(trainLbls),net.Layers,options);
        else
            net=trainNetwork(trainData',categorical(trainLbls),layers,options);
        end
        
        % Generate validation lbls
        valLblsEst=classify(net,valData);
        trainLblsEst=classify(net,trainData);
        valBAcc(currEpoch)=testAcc(valLbls,double(valLblsEst)-1);
        fltrdValBAcc(currEpoch)=median(valBAcc(max(1,currEpoch-4):currEpoch));
        if fltrdValBAcc(currEpoch)==max(fltrdValBAcc)
            bestNet=net;
        end
        trainBAcc=testAcc(trainLbls,double(trainLblsEst)-1);
        fprintf('Epoch: %d/%d, train BAcc: %0.2f, validation BAcc: %0.2f\n',currEpoch,maxEpochs,trainBAcc,fltrdValBAcc(currEpoch));
        
        % Test for ending conditions
        currEpoch=currEpoch+1;
        if ((currEpoch==maxEpochs)||(find(fltrdValBAcc==max(fltrdValBAcc),1,'last')+validationChecks<currEpoch))
            break
        end
    end
    net=bestNet;
    
    % Generate test lbls
    trainLblsEst=classify(net,trainData);
    testLblsEst=classify(net,testData);
    trialLblsEst(C.test(currFold))=double(testLblsEst)-1;
    
    % Compute resulting accuracies on training and testing sets
    trainBAcc=testAcc(trainLbls,double(trainLblsEst)-1);
    testBAcc=testAcc(testLbls,double(testLblsEst)-1);
    fprintf('Train BAcc: %0.2f; validation BAcc: %0.2f; test BAcc: %0.2f\n\n',trainBAcc,max(fltrdValBAcc),testBAcc);
    clear net
end
fprintf('Global test BAcc: %0.2f\n',testAcc(trialLbls,trialLblsEst));

% % Create network
% nNeurons1= 50;
% nNeurons2= 7;
% maxEpochs = 100;
% validationChecks=15;
% miniBatchSize = round(length(trialLbls)/5);
% inputSize = size(trialData{1},1);
% trainErr=zeros(length(nNeurons1),1); %#ok<*PREALL>
% testErr=zeros(length(nNeurons1),1);
% 
% % Define network structure and training options
% layers=[ ...
%     sequenceInputLayer(inputSize)
%     lstmLayer(nNeurons1,'OutputMode','last')
%     fullyConnectedLayer(nNeurons2)
%     dropoutLayer(0.5)
%     fullyConnectedLayer(1)
%     regressionLayer];
% options = trainingOptions('adam', ...
%     'MaxEpochs',1, ...
%     'MiniBatchSize',miniBatchSize, ...
%     'InitialLearnRate',0.01, ...
%     'GradientThreshold',1, ...
%     'Shuffle','never', ...
%     'Verbose',0, ...
%     'ExecutionEnvironment','gpu');
% 
% % Define blocks for k-fold validation
% C.NumTestSets=5;
% % Random grouping
% % C.groups=ceil(rand(size(trialLbls))*C.NumTestSets);
% % Subsequent grouping
% C.groups=ceil(linspace(1/length(trialLbls),C.NumTestSets,length(trialLbls)));
% C.training=@(currGroup)C.groups~=currGroup;
% C.test=@(currGroup)C.groups==currGroup;
% trialLblsEst=zeros(size(trialLbls));
% 
% computeErr=@(x,y)mean((x-y).^2);
% testAcc=@(x,y)(sum((x==1).*(y==1))./sum(x==1)+sum((x==0).*(y==0))./sum(x==0))*.5;
% for currFold=1:C.NumTestSets
%     % Recover training and testing sets
%     trainData=trialData(C.training(currFold));
%     trainLbls=trialLbls(C.training(currFold));
%     testData=trialData(C.test(currFold));
%     testLbls=trialLbls(C.test(currFold));
%     
%     % Split training set in training and validation
%     valSetProportion=.5;
%     rndOrdr=randperm(length(trainData));
%     valIdx=rndOrdr(1:floor(length(trainData)*valSetProportion));
%     trainIdx=setdiff(1:length(trainLbls),valIdx);
%     
%     valData=trainData(valIdx);
%     valLbls=trainLbls(valIdx);
%     trainData=trainData(trainIdx);
%     trainLbls=trainLbls(trainIdx);
%     
%     % Train network
%     currEpoch=1;
%     valErr=zeros(maxEpochs,1);
%     fltrdValErr=ones(maxEpochs,1)*Inf;
%     while true
%         if exist('net','var')
%             net=trainNetwork(trainData',trainLbls,net.Layers,options);
%         else
%             net=trainNetwork(trainData',trainLbls,layers,options);
%         end
%         
%         % Generate validation lbls
%         valLblsEst=predict(net,valData);
%         trainLblsEst=predict(net,trainData);
%         valErr(currEpoch)=computeErr(valLbls,valLblsEst);
%         fltrdValErr(currEpoch)=median(valErr(max(1,currEpoch-4):currEpoch));
%         if fltrdValErr(currEpoch)==min(fltrdValErr)
%             bestNet=net;
%         end
%         trainErr=computeErr(trainLbls,trainLblsEst);
%         fprintf('Epoch: %d/%d, train Err: %0.2f, validation Err: %0.2f\n',currEpoch,maxEpochs,trainErr,min(fltrdValErr));
%         
%         % Test for ending conditions
%         currEpoch=currEpoch+1;
%         if ((currEpoch==maxEpochs)||(find(fltrdValErr==min(fltrdValErr),1,'last')+validationChecks<currEpoch))
%             break
%         end
%     end
%     net=bestNet;
%     
%     % Generate test lbls
%     trainLblsEst=predict(net,trainData);
%     testLblsEst=predict(net,testData);
%     trialLblsEst(C.test(currFold))=testLblsEst;
%     
%     % Compute resulting accuracies on training and testing sets
%     trainErr=computeErr(trainLbls,trainLblsEst);
%     testErr=computeErr(testLbls,testLblsEst);
%     fprintf('Train Err: %0.2f; validation Err: %0.2f; test Err: %0.2f\n\n',trainErr,min(fltrdValErr),testErr);
%     clear net
% end