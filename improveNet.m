close all
clear
clc

% load('C:\Data\2018_05_keyboard\Tests\20180604T161434_longTraining.mat')
% load('fileList.mat');
% obj=MI_MIkeyboard.joinSessions(f);
% load('C:\Data\2018_06_MIkeyboard\S_Chang\20180614T102919.mat')
load('lastNet.mat','net');
f=getFileList;
obj=MI_MIkeyboard.joinSessions(f);

% Recover correct selCounter from UDP transmissions
selCounter=cellfun(@(x)str2double(char(x)),obj.outputLog.keyUDPlogOut);

% Recover beginning and end of key selections from selCounter, as well as
% labels
flexPoints=find((selCounter(1:end-2)-selCounter(2:end-1)).*(selCounter(3:end)-selCounter(2:end-1))>=0)+1;
trialStart=obj.outputLog.time(flexPoints(selCounter(flexPoints)<selCounter(flexPoints+1)));
trialEnd=flexPoints(selCounter(flexPoints)>selCounter(flexPoints-1));
trialLbls=(selCounter(trialEnd)==1)';
trialEnd=obj.outputLog.time(trialEnd);

% Divide collected data in different trials
artifactLimit=15;
timeWindowLength=.2;
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
trialLbls=double(trialLbls);

% Setup parameters
maxEpochs = 100;
validationChecks=15;
miniBatchSize = round(length(trialLbls)/5);
options = trainingOptions('adam', ...
    'ExecutionEnvironment','auto', ...
    'MaxEpochs',1, ...
    'MiniBatchSize',miniBatchSize, ...
    'GradientThreshold',1, ...
    'Verbose',0);

computeBAcc=@(x,y)mean((x-y).^2);
testAcc=@(x,y)(sum((x==1).*(y==1))./sum(x==1)+sum((x==0).*(y==0))./sum(x==0))*.5;
    
% Split training set in training and validation
valSetProportion=.5;
rndOrdr=randperm(length(trialData));
valIdx=rndOrdr(1:floor(length(trialData)*valSetProportion));
trainIdx=setdiff(1:length(trialLbls),valIdx);

valData=trialData(valIdx);
valLbls=trialLbls(valIdx);
trainData=trialData(trainIdx);
trainLbls=trialLbls(trainIdx);

% Train network
currEpoch=1;
valBAcc=zeros(maxEpochs,1);
fltrdValBAcc=zeros(maxEpochs,1);
while true
    net=trainNetwork(trainData',categorical(trainLbls),net.Layers,options);
    
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

save('lastNet.mat','net');
