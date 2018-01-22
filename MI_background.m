classdef MI_background
    properties
        fs=512; % WARNING: DO NOT change this. Amplifier acquisition rate can only be changed from its panel in Simulink
        rawData;
        modelName;
        timeTriggeredEvents;
        outputLog;
        bufferData;
        linearMap;
        figHandle;
        UDPchannels;
        selCounter=0; % Current value of selection counter (selection occurs at 1)
        selThreshold=.8; % Classifier output must be at least this high before selCounter is updated
        selDecay=.1; % Decay rate of selCounter when gaze is outside keys (e.g. .1 means 10% decay per iteration)
        selWait=0; % Training will wait this many seconds after gaze is moved into key before considering moro imagery started
    end
    properties (Dependent)
        currTime;
        actualTarget;
        inTarget;
        currKeyID;
    end
    properties (Hidden)
        CDlength;
        isExpClosed=0;
        isTraining=0;
        selStartTime;
        lastKeyID;
    end
    methods
        %% Constructor
        function obj=MI_background(varargin)
            % If an argument is passed, it must be a structure matching
            % linearMap template (i.e. basically, an output of another
            % session). Some parameters (e.g. sampling frequency of
            % amplifier and buffer window length) cannot be changed
            % directly from here, so please do make sure that they're
            % matching with the current settings of relevant Simulink
            % model.
            
            % Set length of initial countdown, in seconds 
            obj.CDlength=15;
            
            % Define expected inputs and output ports for udp communication
            obj.UDPchannels.outPort=7500;
            obj.UDPchannels.inPort=6100;
                       
            % Normalization buffer length
            obj.bufferData.bufferLength=30; % Data will be normalized in mean and SD over this amount of seconds
            
            % Default parameters
            if nargin==0
                linearMapDefault.fs=obj.fs;
                linearMapDefault.ARmodelOrder=14; % From: "Noninvasive Electroencephalogram Based Control of a Robotic Arm for Reach and Grasp Tasks"
                %             linearMapDefault.bandLims=[8,12,18,25]; % Limits of band of interest - i.e. 8-to-12 and 18-to-25
                %             linearMapDefault.bandLims=ceil((1:.5:24.5)); % Limits of band of interest - single Hertz bands from 1 to 25
                linearMapDefault.bandLims=[10,14]; % Cited work
                % WARNING: the following line DOESN'T DO ANYTHING; it is just a
                % reference. Buffer length has to be changed in the Simulink
                % model (specifically, in the triggeredBuffer mask
                % initialization)
                linearMapDefault.winLength=.4; % Window length, in seconds. Also, cited work
                linearMapDefault.relChannels=1:16; % Use all channels
%                 linearMapDefault.relChannels=[7,11]; % C3 and C4 in 16 el setup
                linearMapDefault.nChannels=length(linearMapDefault.relChannels);
                obj.linearMap=linearMapDefault;
            else
                obj.linearMap=varargin{1};
            end
            
            % Initialize a few things
            obj.outputLog.time=[];
            obj.outputLog.feats=[];
            obj.outputLog.isTraining=[];
            obj.outputLog.isInTarget=[];
            obj.outputLog.selCounter=[];
            obj.outputLog.UDPlog=cell(0);
            obj.outputLog.selStart=[];
            obj.outputLog.selEnd=[];
        end
        
        % Other methods
        function obj=start(obj)
            % Variables on base workspace will be used to trigger closing
            % of experiment
            assignin('base','isExpClosing',0);
            
            % Open udp channels for communication
            obj=openUDPchannels(obj);

            % Opens figure as background. Useful for accepting key presses
            obj=createExpFigure(obj);
            
            % Launch keyboard program
%             !C:\Code\keyboard2\keyboad.exe &
                        
            % Sets name of Simulink model to be used for acquisition
            obj.modelName='SimpleAcquisition_16ch_2014a_RT';
            
            % Prepares Simulink model (i.e. starts recording, basically)
            obj.prepareSimulinkModel;
            
            % Wait for amplifiers to set
            pause(obj.CDlength);
            
            % Generates array of time triggered events
            obj.timeTriggeredEvents{1}=timeTriggeredEvent('expCallback',0);
            obj.timeTriggeredEvents{2}=timeTriggeredEvent('toggleTraining',0);
            
            % Perform bulk of experiment
            obj=manageExperiment(obj);
            
            % Closes exp window and saves data
            obj.closeExp;
        end
        
        function obj=manageExperiment(obj)
            % Generate file name used to save experiment data
            fileName=datestr(now,30);
            while ~evalin('base','isExpClosing')
                pause(0.001);
                for currTTevent=1:length(obj.timeTriggeredEvents);
                    obj=checkAndExecute(obj.timeTriggeredEvents{currTTevent},obj.currTime,obj);
                    pause(0.001);
                end
            end
            obj.isExpClosed=1;
            delete(gcf);
            set_param(obj.modelName,'SimulationCommand','Stop');
            set_param(obj.modelName,'StartFcn','')
            obj.rawData=evalin('base','rawData');
            save(fileName,'obj');
            
            % Release receiver UDP port
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
%             fprintf('%s\n',obj.UDPchannels.lastInput)
            if ~obj.lastKeyID==obj.currKeyID
                obj.selStartTime=obj.currTime;
                obj.selCounter=0;
            end
                        
            % Only execute the following if first EEG data-block has
            % arrived
            if evalin('base','exist(''currData'',''var'')')
                % Recover data buffer from base workspace (Simulink puts them
                % there)
                dataWindow=evalin('base','currData');
                dataTimeStamp=obj.currTime;
                
                % Update normalization buffer and normalize data
                obj=updateBufferData(obj,dataWindow,dataTimeStamp);
                fixDim=@(x)repmat(x,size(dataWindow,1),1);
                if obj.bufferData.SD>0 % Skip first window to prevent Infs
                    dataWindow=(dataWindow-fixDim(obj.bufferData.mean))./fixDim(obj.bufferData.SD);
                end
                
                % If this is first iteration, compute laplacian filter weights
                if ~isfield(obj.linearMap,'lapFltrWeights')
                    [~,obj.linearMap.lapFltrWeights]=MI_background.applyLapFilter(dataWindow);
                end
                
                % Recover bandpower data from data buffer
                BP=MI_background.preprocData(dataWindow,obj.linearMap);
                
                % If training is ongoing, update linear map
                if obj.isTraining
                    if obj.inTarget
                        if obj.currTime-obj.selStartTime>obj.selWait
                            obj=obj.computeLinearMap(BP);
                        end
                    else
                        obj=obj.computeLinearMap(BP);
                    end
                end
                               
                % If current key has already been selected, reset counter to zero
                if obj.selCounter>1 
                    % Communicate whether threshold has been reached
%                     outString=uint8(sprintf('%d',obj.currKeyID));
                    outString=uint8(sprintf('%d',1));

                    obj.selCounter=0;
                    obj.outputLog.selStart=cat(1,obj.outputLog.selStart,obj.selStartTime);
                    obj.outputLog.selEnd=cat(1,obj.outputLog.selEnd,obj.currTime);
                else
                    outString=uint8('0');
                end
                obj.UDPchannels.udps.step(outString);
                fprintf('%s\n',outString)
                
                % If cursor is within target, fills up selection counter
                if obj.inTarget
                    % Use EEG data when available, otherwise use gaze alone
                    if isfield(obj.linearMap,'mat')
                        feats=cat(2,1,BP);
                        currEst=1./(1+exp(-(feats*obj.linearMap.mat)));
%                         fprintf('%0.2f\n',currEst)
                        obj.selCounter=max(0,obj.selCounter+max(0,currEst-obj.selThreshold));
                    else
                        obj.selCounter=obj.selCounter+.1;
                    end
                else
%                     obj.selCounter=obj.selCounter*(1-obj.selDecay);
                end
%                 fprintf('%0.2f\n',obj.selCounter)
                
                % Add relevant info to log
                obj.outputLog.feats=cat(1,obj.outputLog.feats,BP);
                obj.outputLog.isTraining=cat(1,obj.outputLog.isTraining,obj.isTraining);
                obj.outputLog.time=cat(1,obj.outputLog.time,obj.currTime);
                obj.outputLog.isInTarget=cat(1,obj.outputLog.isInTarget,obj.inTarget);
                obj.outputLog.selCounter=cat(1,obj.outputLog.selCounter,obj.selCounter);
                obj.outputLog.UDPlog{end+1}=obj.UDPchannels.lastInput;
            end
            
            % Set next evaluation time for this function
            obj.timeTriggeredEvents{1}.triggersLog=[obj.timeTriggeredEvents{1}.triggersLog,obj.currTime];
            obj.timeTriggeredEvents{1}.nextTrigger=obj.currTime+.05;
        end
        
        function obj=updateBufferData(obj,dataWindow,dataTimeStamp)
            % At the moment, I am only taking a single time point for each
            % dataWindow to estimate mean and sd in the last thirty
            % seconds. This is VERY approximate, better ideas for a
            % different solution are welcome
            if isfield(obj.bufferData,'data')
                obj.bufferData.data=cat(1,obj.bufferData.data,dataWindow(1,:));
                obj.bufferData.timeStamps=cat(1,obj.bufferData.timeStamps,dataTimeStamp);
            else
                obj.bufferData.data=dataWindow(1,:);
                obj.bufferData.timeStamps=dataTimeStamp;
            end
            toBeRemoved=obj.currTime-obj.bufferData.timeStamps>obj.bufferData.bufferLength; % Last value is buffer length, in seconds
            obj.bufferData.data(toBeRemoved,:)=[];
            obj.bufferData.timeStamps(toBeRemoved)=[];
            obj.bufferData.mean=median(obj.bufferData.data,1);
            obj.bufferData.SD=1.4826*median(abs(obj.bufferData.data-repmat(obj.bufferData.mean,size(obj.bufferData.data,1),1)));
%             obj.bufferData.SD=std(obj.bufferData.data,[],1);
        end
        
        function obj=computeLinearMap(obj,BP)
            obj.linearMap.feats{obj.actualTarget}=cat(2,1,BP);
            
            % Define link function
            if ~isfield(obj.linearMap,'mat')
                obj.linearMap.mat=zeros(length(obj.linearMap.feats{obj.actualTarget}'),1);
            end
            
            % Compute output of current classifier
            currEst=1./(1+exp(-(obj.linearMap.feats{obj.actualTarget}*obj.linearMap.mat)));
            
            % Update weights of GLM model
            obj.linearMap.mat=MI_background.updateWeights(obj.linearMap.mat,BP,currEst,obj.inTarget,1e-3);
        end
        
        function obj=toggleTraining(obj)
            if evalin('base','exist(''toggleTraining'',''var'')')&&evalin('base','toggleTraining')
                obj.isTraining=~obj.isTraining;
                assignin('base','toggleTraining',0);
                figure(obj.figHandle)
                if obj.isTraining
                    textHandle=text(-.4,.5,'Training on');
                else
                    textHandle=text(-.4,.5,'Training off');
                end
                set(textHandle,'Color','black','FontSize',28);
            	wait(obj,.5);
                delete(textHandle);
            end
            
            % Set next evaluation time for this function
            obj.timeTriggeredEvents{2}.triggersLog=[obj.timeTriggeredEvents{2}.triggersLog,obj.currTime];
            obj.timeTriggeredEvents{2}.nextTrigger=obj.currTime+.5;
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
        
        function wait(obj,pauseLength)
            startTime=get_param(obj.modelName,'SimulationTime');
            while strcmp(get_param(obj.modelName,'SimulationStatus'),'running')&&get_param(obj.modelName,'SimulationTime')<=startTime+pauseLength
                pause(1/(2*obj.fs));
            end
        end
        
        function obj=openUDPchannels(obj)
            obj.UDPchannels.udpr=dsp.UDPReceiver('LocalIPPort',obj.UDPchannels.inPort,'ReceiveBufferSize',1);
            obj.UDPchannels.udps=dsp.UDPSender('RemoteIPPort',obj.UDPchannels.outPort);
        end
        
        %% Dependent properties
        function cTime=get.currTime(obj)
            if obj.isExpClosed
                cTime=obj.rawData.Time(end);
            else
                cTime=get_param(obj.modelName,'SimulationTime');
            end
        end
        
        function aTarget=get.actualTarget(obj)
                aTarget=double(obj.inTarget)+1;
        end
        
        function inT=get.inTarget(obj) %#ok<MANU>
            inT=1;
        end
        
        function keyID=get.currKeyID(obj)
            if obj.inTarget&&isfield(obj.UDPchannels,'lastInput')&&~isempty(obj.UDPchannels.lastInput)
                keyID=str2double(obj.UDPchannels.lastInput);
            else
                keyID=0;
            end
%             fprintf('%d\n',keyID)
        end
    end
    methods (Static)
        function preprocessedData=preprocData(dataWindow,linearMap)
            preprocessedData=zeros(1,linearMap.nChannels*length(linearMap.bandLims)/2); % Initialize band power matrix.
            lapData=dataWindow*linearMap.lapFltrWeights;
            for currCh=1:linearMap.nChannels
                % Compute power in input window using Yule-Walker PSD
                pxx=pyulear(detrend(lapData(:,currCh)),linearMap.ARmodelOrder);
                
                % Compute power in bands of interest
                binCenters=linspace(1/linearMap.fs,linearMap.fs/2,size(pxx,1));
                for currBand=1:length(linearMap.bandLims)/2
                    [~,bandStart]=min(abs(binCenters-linearMap.bandLims(currBand*2-1)));
                    [~,bandEnd]=min(abs(binCenters-linearMap.bandLims(currBand*2)));
                    preprocessedData(:,(currCh-1)*length(linearMap.bandLims)/2+currBand)=sum(pxx(bandStart:bandEnd,:))';
                end
            end
            preprocessedData=log(preprocessedData);
        end
        
        function [outData,fltrWeights]=applyLapFilter(inData)
            try
                load('elMap16.mat')
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
                
        function wOut=updateWeights(wIn,feats,E,t,lr)
            % feats new sample
            % E current classifier prediction
            % t true label
            feats=[1,feats];
            wOut=wIn+((lr*(t-E))'*feats)';
        end
        
        function [coords,inKey,keyID]=extractLogs(log)
            coords=zeros(length(log),3);
            inKey=zeros(length(log),1);
            keyID=cell(length(log),1);
            for currLog=1:length(log)
                if isempty(log{currLog})
                    coords(currLog,:)=0;
                    inKey(currLog,:)=0;
                    keyID='';
                else
                    cBracePos=strfind(log{currLog},'}');
                    coords(currLog,:)=str2num(log{currLog}(2:cBracePos-1)); %#ok<ST2NM>
                    inKey(currLog)=numel(strfind(log{currLog},'True'))>0;
                    semiColPos=strfind(log{currLog},';');
                    keyID{currLog}=log{currLog}(semiColPos(end)+1:end);
                end
            end
        end
    end
end

function simulinkModelStartFcn(modelName) %#ok<DEFNU>
% Start function for Simulink model.
blockName=sprintf('%s/triggeredBuffer/Buffer',modelName);
assignin('base','listener',add_exec_event_listener(blockName,'PostOutputs',@acquireBufferedData));
end

function acquireBufferedData(block,~)
assignin('base','currData',block.OutputPort(1).Data);
assignin('base','currTime',block.SampleTime);
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
    MI_background.closeExp;
end
if strcmp(eventdata.Key,'p')
    keyboard;
    %     assignin('base','pauseNextTrial',1)
end
if strcmp(eventdata.Key,'t')
    assignin('base','toggleTraining',1);
end
end

function OnClosing(~,~)
% Overrides normal closing procedure so that regardless of how figure is
% closed logged data is not lost
MI_background.closeExp;
end