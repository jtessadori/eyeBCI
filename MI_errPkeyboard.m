classdef MI_errPkeyboard
    properties
        fs=512; % WARNING: DO NOT change this. Amplifier acquisition rate can only be changed from its panel in Simulink, it is here for reference
        errPbufferLength=5; % Same as above. Parameters can be changed in triggeredBuffer block initialization mask in the Simulink block
        rawData;
        condition;
        modelName;
        timeTriggeredEvents;
        errPclassifier;
        outputLog;
        timingParams;
        UDPchannels;
        lastKeyID;
        selStartTime;
        selCounter=0;
        selDuration=1.5;
        selSound;
        keyPressTimeBuffer=[];
        errorInterval=10; % Force an error approximately once every this many seconds during training
        lastErrorTime=30;
        targetErrorRate=0.05;
        estAvg=0.05;
    end
    properties (Dependent)
        currTime;
        inTarget;
        currKeyID;
        ETcursorPos;
    end
    properties (Hidden)
        CDlength;
        isExpClosed=0;
        isTraining;
        errPdataWindow;
        isPausing=0;
        figHandle;
        possibleConditions={'Training','Test'};
        selDurationBak;
    end
    methods
        %% Constructor
        function obj=MI_errPkeyboard(varargin)
            % Exactly one argument may be passed: another instance of a
            % MI_errPkeyboard class. The only properties that will be used are the
            % classifiers, everything else is left at default values.
            if nargin>=1
                obj.errPclassifier=varargin{1}.errPclassifier;
            end
            % Some parameters (e.g. sampling frequency of amplifier and
            % buffer window length) cannot be changed directly from here,
            % so please do make sure that they're matching with the current
            % settings of relevant Simulink model.
            
            % Set length of initial countdown, in seconds
            obj.CDlength=15;
            
            % Define expected inputs and output ports for udp communication
            obj.UDPchannels.outPort=9000;
            obj.UDPchannels.inPort=8051;
            
            % Define timing parameters
            obj.timingParams.errPestimationLength=1; % Length of window used for errP estimation, following movement
            
            % Define sound to be played at each letter selection
            obj.selSound.f=32000;
            t=linspace(0,1,obj.selSound.f);
            obj.selSound.y=betapdf(t,2,20).*sin(2*pi*t*700);
                         
            % Initialize log entries
            % Each time step
            obj.outputLog.time=[];
            obj.outputLog.cursorPos=[];
            obj.outputLog.isTraining=[];
            obj.outputLog.isInTarget=[];
            obj.outputLog.selCounter=[];
            obj.outputLog.UDPlogIn=[];
            obj.outputLog.UDPlogOut=[];
            % Following key selection
            obj.outputLog.selStart=[];
            obj.outputLog.selEnd=[];
            obj.outputLog.selID=[];
            obj.outputLog.errPtimes=[];
            obj.outputLog.errPest=[];
            obj.outputLog.errPfeats=[];
            obj.outputLog.forcedError=[];
            
            % Ask user whether to start experiment right away
            clc;
            if ~strcmpi(input('Start experiment now? [Y/n]\n','s'),'n')
                obj=runExperiment(obj);
            end
        end
        
        % Other methods
        function obj=runExperiment(obj)
            % Variables on base workspace will be used to trigger closing
            % of experiment
            assignin('base','isExpClosing',0);
            
            % Sets name of Simulink model to be used for acquisition
            obj.modelName='SimpleAcquisition_16ch_2014a_RT_preProc';
            
            % Prompts user to select a condition
            obj=selectCondition(obj);
            obj=setConditionSpecificParams(obj);
            
            % Open udp channels for communication
            obj=openUDPchannels(obj);

            % Opens figure as background. Useful for accepting key presses
            obj=createExpFigure(obj);
            
            % Launch keyboard program
%             !C:\Code\keyboard2\keyboad.exe &
            !C:\Code\keyboard_jacopo_new\keyboad.exe &
                                    
            % Prepares Simulink model (i.e. starts recording, basically)
            obj.prepareSimulinkModel;
            
            % Generates array of time triggered events
            obj.timeTriggeredEvents{1}=timeTriggeredEvent('expCallback',obj.CDlength);
            obj.timeTriggeredEvents{2}=timeTriggeredEvent('processErrPbuffer',Inf);
            
            % Shows a countdown
            obj.startCountdown(obj.CDlength);
            
            % Perform bulk of experiment
            obj=manageExperiment(obj);
            
            % Closes exp window and saves data
            obj.closeExp;
        end
        
        function obj=manageExperiment(obj)
            % Generate file name used to save experiment data
            fileName=datestr(now,30);
                        
            % Set pausing variable on base workspace
            assignin('base','togglePause',0);
            
            % Experiment control loop
            while ~evalin('base','isExpClosing')
                pause(0.001);
                for currTTevent=1:length(obj.timeTriggeredEvents)
                    obj=checkAndExecute(obj.timeTriggeredEvents{currTTevent},obj.currTime,obj);
                    pause(0.0001);
                end
                if (evalin('base','togglePause'))
                    assignin('base','togglePause',0);
                    if obj.isPausing
                        obj.isPausing=0;
                        set_param(obj.modelName,'SimulationCommand','Continue');
                    else
                        obj.isPausing=1;
                        set_param(obj.modelName,'SimulationCommand','Pause');
                    end
                end
            end
            pause(1);
            obj.isExpClosed=1;
            delete(gcf);
            set_param(obj.modelName,'SimulationCommand','Stop');
            set_param(obj.modelName,'StartFcn','')
            obj.rawData=evalin('base','rawData');
            save(fileName,'obj');
            
            % Release receiver UDP port
            obj.UDPchannels.udpr.release;
            obj.UDPchannels.udpr.isDone;
            
            % Clear variables from base workspace
            evalin('base','clear listener*');
            evalin('base','clear toggleTraining');
        end
        
        function obj=createExpFigure(obj)
            % Set figure properties
            obj.figHandle=figure;
            set(obj.figHandle,'Tag',mfilename,...
                'Toolbar','none',...
                'MenuBar','none',...
                'Units','normalized',...
                'Resize','off',...
                'NumberTitle','off',...
                'Name','',...
                'WindowKeyPressFcn',@KeyPressed,...
                'CloseRequestFcn',@OnClosing,...
                'WindowButtonMotionFcn',@onMouseMove);
            
            % Remove figure axis
            axis([-1,1,-1,1]);
            axis('off')
        end
        
        function obj=expCallback(obj)
            % Recover incoming data and update relevant properties
            obj.lastKeyID=obj.currKeyID;
            obj.UDPchannels.lastInput=char(obj.UDPchannels.udpr.step)';
            if ~strcmp(obj.lastKeyID,obj.currKeyID)
                obj.selStartTime=obj.currTime;
                obj.selCounter=0;
            end
            
            % If cursor is within target, fills up selection counter
            if obj.inTarget
                obj.selCounter=(obj.currTime-obj.selStartTime)/obj.selDuration;
            else
                obj.selCounter=0;
            end
                                                                                                                  
            % If current key has already been selected, reset counter to zero
            if obj.selCounter>=1
                % Play sound
                sound(obj.selSound.y,obj.selSound.f);
                
                % Communicate whether threshold has been reached
                if obj.isTraining && rand<.1%obj.currTime>obj.lastErrorTime+obj.errorInterval
                    % If applicable, commit a mistake
                    outString=uint8(sprintf('%d',3));
                    obj.outputLog.forcedError=cat(1,obj.outputLog.forcedError,1);
                    obj.lastErrorTime=obj.currTime;
                else
                    outString=uint8(sprintf('%d',1));
                    obj.outputLog.forcedError=cat(1,obj.outputLog.forcedError,0);
                end
                
                % Reset selection counter and update relevant logs
                obj.selCounter=0;
                obj.outputLog.selStart=cat(1,obj.outputLog.selStart,obj.selStartTime);
                obj.outputLog.selEnd=cat(1,obj.outputLog.selEnd,obj.currTime);
                obj.outputLog.selID{end+1}=obj.currKeyID;
                obj.selStartTime=obj.currTime;
                
                % Update buffer of key presses
                obj.keyPressTimeBuffer=[obj.keyPressTimeBuffer,obj.currTime];
            else                
                outString=uint8(sprintf('%0.2f',min(obj.selCounter,0.989))); % Threshold output to prevent sending 1 as result of a rounding error
            end
            obj.UDPchannels.udps.step(outString);
            
            % Add relevant info to log
            obj.outputLog.time=cat(1,obj.outputLog.time,obj.currTime);
            obj.outputLog.isTraining=cat(1,obj.outputLog.isTraining,obj.isTraining);
            obj.outputLog.isInTarget=cat(1,obj.outputLog.isInTarget,obj.inTarget);
            obj.outputLog.selCounter=cat(1,obj.outputLog.selCounter,obj.selCounter);
            obj.outputLog.UDPlogIn{end+1}=obj.UDPchannels.lastInput;
            obj.outputLog.UDPlogOut{end+1}=outString;
            obj.outputLog.cursorPos=cat(1,obj.outputLog.cursorPos,obj.ETcursorPos);
            
            % Set errP evaluation 1s from now
            if ~isempty(obj.keyPressTimeBuffer)
                obj.timeTriggeredEvents{2}.nextTrigger=obj.keyPressTimeBuffer(1)+1;
            end
            
            % Set next evaluation time for this function
            obj.timeTriggeredEvents{1}.triggersLog=[obj.timeTriggeredEvents{1}.triggersLog,obj.currTime];
            obj.timeTriggeredEvents{1}.nextTrigger=obj.currTime;
        end
                
        function obj=processErrPbuffer(obj)
            % Log timing and reset trigger time
            obj.timeTriggeredEvents{2}.nextTrigger=Inf;
            obj.timeTriggeredEvents{2}.triggersLog=[obj.timeTriggeredEvents{2}.triggersLog,obj.currTime];
            
            % Recover data buffer from base workspace (Simulink puts them
            % there)
            dataWindow=evalin('base','currErrPdata');
            dataTimeStamp=obj.currTime;
            
            % Data buffer is longer than required and not exactly synched
            % with cursor movement. Recover correct section (starting from
            % key selection)
            relKeySelection=obj.keyPressTimeBuffer(1);
            obj.keyPressTimeBuffer(1)=[];
            currDelay=dataTimeStamp-relKeySelection;
            dataWindow=dataWindow(round((obj.errPbufferLength-currDelay)*obj.fs+1):round((obj.errPbufferLength-currDelay+obj.timingParams.errPestimationLength)*obj.fs),:);
            
            % Recover frequency features
            freqFeats=evalin('base','BP');
            
            % Join freq and time Data
            feats=cat(1,freqFeats,reshape(resample(dataWindow,64,512),[],1));
            
            % If a classifier is available, perform classification
            if ~isempty(obj.errPclassifier)
                currEst=obj.errPclassifier.clsfr.predict(feats(obj.errPclassifier.featsIdx)');
            else
                currEst=NaN;
            end
            
            % If prediction is finite and not in training, update sel
            % duration
            if ~obj.isTraining && isfinite(currEst)
                obj.estAvg=0.99*obj.estAvg+0.01*currEst;
                obj.selDuration=max(obj.selDuration*(1-obj.estAvg+obj.targetErrorRate),0.2);
            end

            % If not training, perform classification
%             if ~obj.isTraining
%                 currEst=obj.errPclassifier.clsfr.predict(feats(obj.errPclassifier.featsIdx)');
%             else
%                 currEst=NaN;
%             end
%             if ~obj.isTraining
%                currEst=1;
%                if currEst>.5 % i.e. an error is detected
%                     obj.selDuration=(obj.selDuration+(obj.durationMin*obj.durationDecay)).*(1-obj.durationDecay);
%                 else
%                     obj.selDuration=obj.selDuration.*(1+obj.durationIncrease);
%                 end
%             else
%                 currEst=NaN;
%             end
            
            % Log relevant data
            obj.outputLog.errPtimes=cat(1,obj.outputLog.errPtimes,dataTimeStamp);
            obj.outputLog.errPest=cat(1,obj.outputLog.errPest,currEst);
            obj.outputLog.errPfeats=cat(1,obj.outputLog.errPfeats,feats');
        end
                                                                
        function obj=selectCondition(obj)
            currCond=0;
            while true
                clc;
                for currPossibleCond=1:length(obj.possibleConditions)
                    fprintf('[%d] - %s;\n',currPossibleCond,obj.possibleConditions{currPossibleCond});
                end
                currCond=input('\nPlease select desired condition: ');
                if ismember(currCond,1:length(obj.possibleConditions))
                    break
                end
            end
            obj.condition.conditionID=currCond;
        end
        
        function obj=setConditionSpecificParams(obj)
            % 'Training' 'Testing'
            switch obj.condition.conditionID
                case 1
                    obj.isTraining=1;
                    obj.selDurationBak=obj.selDuration;
                case 2
                    obj.isTraining=0;
            end
        end
        
        function prepareSimulinkModel(obj)
            % Check whether simulink model file can be found
            if ~exist(obj.modelName,'file')
                warning('Cannot find model %s.\nPress Enter to continue.\n',obj.modelName);
                input('');
                [fileName,pathName]=uigetfile('*.slx','Select Simulink model to load:');
                obj.modelName=sprintf('%s\\%s',pathName,fileName);
            end
            % Load model
            load_system(obj.modelName);
            
            % Check whether simulation was already running, and, in case,
            % stop it
            if bdIsLoaded(obj.modelName)&&strcmp(get_param(obj.modelName,'SimulationStatus'),'running')
                set_param(obj.modelName,'SimulationCommand','Stop');
            end
            
            % Add event listener to triggered buffer event.
            set_param(obj.modelName,'StartFcn',sprintf('simulinkModelStartFcn(''%s'')',obj.modelName))
            set_param(obj.modelName,'StopTime','inf');
            set_param(obj.modelName,'FixedStep',['1/',num2str(obj.fs)]);
            set_param(obj.modelName,'SimulationCommand','Start');
        end
        
        function obj=openUDPchannels(obj)
            obj.UDPchannels.udpr=dsp.UDPReceiver('LocalIPPort',obj.UDPchannels.inPort,'ReceiveBufferSize',1);
            obj.UDPchannels.udps=dsp.UDPSender('RemoteIPPort',obj.UDPchannels.outPort);
        end
                
        function wait(obj,pauseLength)
            startTime=get_param(obj.modelName,'SimulationTime');
            while strcmp(get_param(obj.modelName,'SimulationStatus'),'running')&&get_param(obj.modelName,'SimulationTime')<=startTime+pauseLength
                pause(1/(2*obj.fs));
            end
        end
        
        function startCountdown(obj,nSecs)
            % countdown to experiment start
            figure(obj.figHandle)
            for cntDown=nSecs:-1:1
                if ~exist('textHandle','var')
                    textHandle=text(-.05,.5,num2str(cntDown));
                else
                    set(textHandle,'String',num2str(cntDown));
                end
                set(textHandle,'Color','white','FontSize',64);
                pause(1);
            end
            delete(textHandle);
        end
        
        function [clsfr,cvAcc]=computeErrPclassifier(obj)
            % Recover feats and labels
            [allFeats,lbls]=recoverErrPdata(obj);
                        
%             % Make a first selection of relevant features
%             classLbls=unique(lbls);
%             m=zeros(length(classLbls),size(allFeats,2));
%             md=zeros(size(m));
%             for currClass=1:length(classLbls)
%                 % Use median and mad as proxy for mean and sd, to reduce
%                 % relevance of artifacts
%                 m(currClass,:)=median(allFeats(lbls==classLbls(currClass),:));
%                 md(currClass,:)=1.4826*mad(allFeats(lbls==classLbls(currClass),:),1);
%             end
%             computeWorth=@(m1,m2,md1,md2)abs((m1-m2)./sqrt(md1.^2+md2.^2));
%             featWorth=computeWorth(m(1,:),m(2,:),md(1,:),md(2,:));
%             
%             % Keep features with a worth greater than 0.3 (keep at least
%             % 15)
%             [sortedWorth,featOrdr]=sort(featWorth,'descend');
%             goodFeatsNumber=sum(sortedWorth>.3);
%             goodFeatsIdx=featOrdr(1:max(15,goodFeatsNumber));
%             feats=allFeats(:,goodFeatsIdx);
            
            % Train classifier and present cross-validation results
            fprintf('Training ErrP classifier. Please be patient, it will take some time...\n\n');
            [clsfr,~,cvAcc]=testClassifier2(lbls,allFeats,'blocktype','subsequent','nblocks',10,'threshold',.4,'selectionType','zScore');
%             [clsfr,~,cvAcc]=testClassifier2(lbls,allFeats,'blocktype','subsequent','nblocks',10,'threshold',.65,'selectionType','histOverlap');
        end
        
        function cvAcc=testErrPrequiredLength(obj)
            % Recover feats and labels
            [allFeats,lbls]=recoverErrPdata(obj);
            
            % Preserve original data
            allFeatsBak=allFeats;
            lblsBak=lbls;
            
            % Test classifier accuracy for different lengths
            nTrials=25:25:length(lblsBak);
            cvAcc=zeros(length(nTrials),1);
            for currLength=1:length(nTrials)
                allFeats=allFeatsBak(1:nTrials(currLength),:);
                lbls=lblsBak(1:nTrials(currLength));
                                
                % Present cross-validation results
%                 [~,~,cvAcc(currLength)]=testClassifier2(lbls,allFeats,'blocktype','subsequent','nblocks',6,'threshold',.2);
                [~,~,cvAcc(currLength)]=testClassifier2(lbls,allFeats,'blocktype','subsequent','nblocks',10,'featureN',4,'selectionType','histOverlap');
                plot(nTrials(1:currLength),cvAcc(1:currLength))
                pause(0.01);
                fprintf('%d/%d\n',currLength,length(nTrials));
            end
        end
        
        function [allFeats,lbls]=recoverErrPdata(obj)
            lbls=obj.outputLog.forcedError;
            allFeats=obj.outputLog.errPfeats;
        end
        
        function BACC=testErrPclassifier(obj)
            % Compute predictions
            errPest=obj.errPclassifier.clsfr.predict(obj.outputLog.errPfeats(:,obj.errPclassifier.featsIdx));
            
            % Compute BACC
            testAcc=@(x,y)(sum((x==1).*(y==1))./sum(x==1)+sum((x==0).*(y==0))./sum(x==0))*.5;
            BACC=testAcc(obj.outputLog.forcedError,errPest);
%             [~,~,~,AUC]=perfcurve(obj.outputLog.correctMovement(1:end-1),errPest,1);
        end
                
        function plotErrPs(obj)
            % Normalize data
            normalize=@(x)(x-repmat(mean(x),size(x,1),1))./repmat(1.4826*mad(x),size(x,1),1);
            normData=normalize(obj.rawData.data);
            
            [B,A]=cheby1(4,6,[1,10]/(obj.fs/2));
            lbls=obj.outputLog.correctMovement;
            
            % Apply spatial and freq filters
%             lapData=obj.applyLapFilter(obj.rawData.data);
            carData=MI_errPkeyboard.applyLapFilter(normData);
            freqData=filter(B,A,carData);
            
            relWins=zeros(length(obj.outputLog.errPtimes),obj.fs*2,size(obj.rawData.data,2));
            for currWin=1:size(relWins,1)
                relWins(currWin,:,:)=freqData((obj.outputLog.errPtimes(currWin)-.5)*obj.fs+1:(obj.outputLog.errPtimes(currWin)+1.5)*obj.fs,:);
            end
            lbls=lbls(1:length(obj.outputLog.errPtimes));
            
            t=linspace(-0.5,1.5,obj.fs*2);
            load('elMap16_MI_err.mat') %#ok<LOAD>
            for currCh=1:16
                subplot(4,4,currCh);
                plot(t,squeeze(median(relWins(lbls==0,:,currCh))),'k');
                hold on;
                plot(t,squeeze(median(relWins(lbls==1,:,currCh))),'r');
                plot(t,squeeze(median(relWins(lbls==0,:,currCh)))-squeeze(median(relWins(lbls==1,:,currCh))),'g','LineWidth',2);
                axis([-.5,1.5,-.1,.1]);
                set(gca,'XTickLabel',[],'YTickLabel',[]);
                xlabel(elMap16.elName{currCh});
            end
        end
                
        function outData=preProcData(obj)
            % Prepare freq filters
            windowLength=.4;
            nBands=19;
            nChannels=16;
            for currFreq=1:nBands
                [B(currFreq,:),A(currFreq,:)]=cheby1(2,6,([2*currFreq,2*(currFreq+1)]/obj.fs)/2); %#ok<AGROW>
            end
            Bfir=ones(1,round(obj.fs*windowLength))/round(obj.fs*windowLength);
            
            % Prepare spatial filters
            fltrWeights=zeros(16);
            try
                load('elMap16_MI_err.mat') %#ok<LOAD>
            catch ME %#ok<NASGU>
                warning('''elMap16_MI_err'' not found. Electrode map required for laplacian filters.');
                return;
            end
            for currEl=1:16
                neighborsMap=zeros(size(elMap16.elMat));
                neighborsMap(elMap16.elMat==currEl)=1;
                neighborsMap=imdilate(neighborsMap,strel('diamond',1));
                neighborsMap(elMap16.elMat==currEl)=0;
                validNeighbors=logical(neighborsMap.*elMap16.elMat);
                fltrWeights(currEl,elMap16.elMat(validNeighbors))=-1/sum(sum(validNeighbors));
                fltrWeights(currEl,currEl)=1;
            end
            
            % Apply spatial filters
            lapData=obj.rawData.data*fltrWeights;
            
            % Apply freq filters
            outData=repmat(obj.rawData.data,1,nBands);
            for currFreq=1:nBands
                outData(:,(currFreq-1)*nChannels+1:currFreq*nChannels)=filter(B(currFreq,:),A(currFreq,:),lapData);
            end
            outData=outData.^2;
            for currFreq=1:nBands
                outData(:,(currFreq-1)*nChannels+1:currFreq*nChannels)=filter(Bfir,1,outData(:,(currFreq-1)*nChannels+1:currFreq*nChannels));
            end
            outData=log10(outData);
        end
        
        function obj=attachPhase(obj,otherSession)
            timeStep=median(diff(obj.rawData.time));
            otherSession.rawData.time=linspace(obj.rawData.time(end)+timeStep,obj.rawData.time(end)+otherSession.rawData.time(end),length(otherSession.rawData.time));
            joiningFields={'cursorPos','errPest','errPfeats','selCounter','isTraining','isInTarget','forcedError'};
            updatingFields={'time','errPtimes','selStart','selEnd'};
            for currUpField=1:length(updatingFields)
                obj.outputLog.(updatingFields{currUpField})=cat(1,obj.outputLog.(updatingFields{currUpField}),otherSession.outputLog.(updatingFields{currUpField})+obj.rawData.time(end));
            end
            for currJoinField=1:length(joiningFields)
                try
                obj.outputLog.(joiningFields{currJoinField})=cat(1,obj.outputLog.(joiningFields{currJoinField}),otherSession.outputLog.(joiningFields{currJoinField}));
                catch
                    keyboard;
                end
            end
            for currTTE=1:length(obj.timeTriggeredEvents)
                obj.timeTriggeredEvents{currTTE}.triggersLog=cat(2,obj.timeTriggeredEvents{currTTE}.triggersLog,otherSession.timeTriggeredEvents{currTTE}.triggersLog+obj.rawData.time(end));
            end
            obj.rawData=append(obj.rawData,otherSession.rawData);
        end
                
        function plotSpectrum(obj)
            % Around 5000 samples, in current setup, are affected with
            % startup artifact
            data=obj.rawData.data(5001:end,:);
            
            f=linspace(1/obj.fs,obj.fs,length(data));
            psd=abs(fft(detrend(data)));
            loglog(f,medfilt1(max(psd,[],2),7));
            xlim([f(2),obj.fs/2]);
        end
        
        %% Dependent properties
        function cTime=get.currTime(obj)
            if obj.isExpClosed
                cTime=obj.rawData.Time(end);
            else
                cTime=get_param(obj.modelName,'SimulationTime');
            end
        end
                
        function inT=get.inTarget(obj)
            inT=0;
            if isfield(obj.UDPchannels,'lastInput')&&~isempty(obj.UDPchannels.lastInput)&&~sum(isnan(obj.UDPchannels.lastInput))
                if numel(strfind(obj.UDPchannels.lastInput,'True'))
                    inT=1;
                end
            end
        end
        
        function keyID=get.currKeyID(obj)
            if obj.inTarget&&isfield(obj.UDPchannels,'lastInput')&&~isempty(obj.UDPchannels.lastInput)&&~sum(isnan(obj.UDPchannels.lastInput))
                semiColumnPos=strfind(obj.UDPchannels.lastInput,';');
                keyID=obj.UDPchannels.lastInput(semiColumnPos(end)+1:end);
            else
                keyID=0;
            end
%             fprintf('%d\n',keyID)
        end
        
        function pos=get.ETcursorPos(obj)
            if isfield(obj.UDPchannels,'lastInput')&&~isempty(obj.UDPchannels.lastInput)&&~sum(isnan(obj.UDPchannels.lastInput))
                lastEntry=obj.UDPchannels.lastInput;
                semiColumnPos=strfind(lastEntry,';');
                lastEntry(semiColumnPos)=',';
                pos=str2num(lastEntry(2:semiColumnPos(2)-1)); %#ok<ST2NM>
            else
                pos=nan(1,2);
            end
        end
    end
    methods (Static)
        function [freqFeats,timeFeats]=preprocessData(dataWins)
            % This function takes either one time window as input (during
            % testing) or a vector of them (during training). Reshape
            % single window to make it consistent
            if length(size(dataWins))==2
                dataWins=reshape(dataWins,1,size(dataWins,1),size(dataWins,2));
            end
            [nWins,~,nChannels]=size(dataWins);
            timeFeats=zeros(size(dataWins,1),round(size(dataWins,2)/8),size(dataWins,3));
            freqFeats=zeros(nWins,129,nChannels);
            % Preprocess each input window
            for currWin=1:nWins
                for currCh=1:nChannels
                    relData=squeeze(dataWins(currWin,:,currCh));
                    % Normalize: set first sample to zero, sd to 1
                    relData=(relData-relData(1))/std(relData);
                    % Remove linear trend
                    relData=detrend(relData);
                    timeFeats(currWin,:,currCh)=resample(relData,64,512); % Resample time features at 64Hz (assuming a 512Hz original sampling rate)
                    % Compute log of bandpower
                    freqFeats(currWin,:,currCh)=pyulear(relData.*blackman(length(relData))',16);
                end                
            end
            % Consider only frequencies up to ~60Hz
            freqFeats(:,31:end,:)=[];
%             % Normalize, then extract logs
%             freqFeats=freqFeats./repmat(sum(freqFeats,3),1,1,size(freqFeats,3));
%             freqFeats=log(freqFeats);
        end
        
        function [outData,fltrWeights]=applyLapFilter(inData)
            try
                load('elMap16_MI_err.mat') %#ok<LOAD>
            catch ME %#ok<NASGU>
                warning('''elMap.mat'' not found. Electrode map required for laplacian filters.');
                outData=[];
                return;
            end
            fltrWeights=zeros(size(inData,2));
            for currEl=1:size(inData,2)
                neighborsMap=zeros(size(elMap16.elMat));
                neighborsMap(elMap16.elMat==currEl)=1;
                neighborsMap=imdilate(neighborsMap,strel('diamond',1));
                neighborsMap(elMap16.elMat==currEl)=0;
                validNeighbors=logical(neighborsMap.*elMap16.elMat);
                fltrWeights(currEl,elMap16.elMat(validNeighbors))=-1/sum(sum(validNeighbors));
                fltrWeights(currEl,currEl)=1;
            end
            outData=inData*fltrWeights';
        end
        
        function closeExp
            % Signals experiment to close
            assignin('base','isExpClosing',1);
        end
                        
        function objLong=joinSessions(fileNames)
            load(fileNames{1});
            objLong=obj;
            for currFile=2:length(fileNames)
                load(fileNames{currFile});
                objLong=attachPhase(objLong,obj);
            end
        end
    end
end

function simulinkModelStartFcn(modelName) %#ok<DEFNU>
% Start function for Simulink model.
blockName=sprintf('%s/filterBlock/log',modelName);
assignin('base','listenerMI',add_exec_event_listener(blockName,'PostOutputs',@acquireFreqFeats));
blockName=sprintf('%s/filterBlock/errP_buffer',modelName);
assignin('base','listenerErrP',add_exec_event_listener(blockName,'PostOutputs',@acquireErrPbufferedData));
end

function acquireFreqFeats(block,~)
assignin('base','BP',block.OutputPort(1).Data);
assignin('base','currTime',block.SampleTime);
end

function acquireErrPbufferedData(block,~)
assignin('base','currErrPdata',block.OutputPort(1).Data);
assignin('base','currErrPtime',block.SampleTime);
end

function onMouseMove(~,~)
% Makes mouse pointer invisible
if ~strcmp(get(gcf,'Pointer'),'custom')
    set(gcf,'PointerShapeCData',NaN(16));
    set(gcf,'Pointer','custom');
end
end

function KeyPressed(~,eventdata,~)
% This is called each time a keyboard key is pressed while the mouse cursor
% is within the window figure area
if strcmp(eventdata.Key,'escape')
    MI_errPkeyboard.closeExp;
end
if strcmp(eventdata.Key,'p')
    keyboard;
    %     assignin('base','pauseNextTrial',1)
end
if strcmp(eventdata.Key,'t')
    assignin('base','toggleTraining',1);
end
if strcmp(eventdata.Key,'z')
    assignin('base','togglePause',1);
end
end

function OnClosing(~,~)
% Overrides normal closing procedure so that regardless of how figure is
% closed logged data is not lost
MI_errPkeyboard.closeExp;
end