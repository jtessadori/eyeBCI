close all
clear
clc

% load('C:\Data\2018_05_keyboard\Tests\20180604T161434_longTraining.mat')
load('C:\Data\2018_06_MIkeyboard\20180608T144422_oldNetTest.mat')

% Recover correct selCounter from UDP transmissions
selCounter=cellfun(@(x)str2double(char(x)),obj.outputLog.keyUDPlogOut);

% Recover beginning and end of key selections from selCounter, as well as
% labels
flexPoints=find((selCounter(1:end-2)-selCounter(2:end-1)).*(selCounter(3:end)-selCounter(2:end-1))>=0)+1;
% flexPoints=union(flexPoints,find(diff(diff(selCounter)==0)==1));
trialStart=obj.outputLog.time(flexPoints(selCounter(flexPoints)<selCounter(flexPoints+1)));
trialEnd=flexPoints(selCounter(flexPoints)>selCounter(flexPoints-1));
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
    relIdx=(obj.rawData.time>relLims(1))&(obj.rawData.time<=relLims(2));
    relData=obj.rawData.data(relIdx,:)';
    trialData{currTrial}=relData;
end
toBeRemoved=cellfun(@(x)isempty(x),trialData);
trialData(toBeRemoved)=[];
trialLbls(toBeRemoved)=[];

% Create network
nNeurons= 100;
maxEpochs = 100;
validationChecks=15;
miniBatchSize = round(length(trialLbls)/5);
inputSize = size(trialData{1},1);
numClasses = 2;

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
    'Verbose',0);

testAcc=@(x,y)(sum((x==1).*(y==1))./sum(x==1)+sum((x==0).*(y==0))./sum(x==0))*.5;
    
% Split training set in training and validation
valSetProportion=.5;
rndOrdr=randperm(length(trialData));
valIdx=rndOrdr(1:floor(length(trialData)*valSetProportion));
trainIdx=setdiff(1:length(trialLbls),valIdx);

valData=trialData(valIdx);
valLbls=trialLbls(valIdx);
trialData=trialData(trainIdx);
trialLbls=trialLbls(trainIdx);

% Train network
currEpoch=1;
valBACC=zeros(maxEpochs,1);
fltrdValBACC=zeros(maxEpochs,1);
while true
    if exist('net','var')
        net=trainNetwork(trialData',categorical(trialLbls),net.Layers,options);
    else
        net=trainNetwork(trialData',categorical(trialLbls),layers,options);
    end
    
    % Generate validation lbls
    valLblsEst=classify(net,valData);
    trialLblsEst=classify(net,trialData);
    valBACC(currEpoch)=testAcc(valLbls,double(valLblsEst)-1);
    fltrdValBACC(currEpoch)=median(valBACC(max(1,currEpoch-4):currEpoch));
    if fltrdValBACC(currEpoch)==max(fltrdValBACC)
        bestNet=net;
    end
    trainBACC=testAcc(trialLbls,double(trialLblsEst)-1);
    fprintf('Epoch: %d/%d, train BACC: %0.2f, validation BACC: %0.2f\n',currEpoch,maxEpochs,trainBACC,fltrdValBACC(currEpoch));
    
    % Test for ending conditions
    currEpoch=currEpoch+1;
    if ((currEpoch==maxEpochs)||(find(fltrdValBACC==max(fltrdValBACC),1,'last')+validationChecks<currEpoch))
        break
    end
end
net=bestNet;