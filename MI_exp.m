classdef MI_exp
    properties
        fs=512; % WARNING: DO NOT change this. Amplifier acquisition rate can only be changed from its panel in Simulink
        rawData;
        modelName;
        timeTriggeredEvents;
        outputLog;
        bufferData;
        linearMap;
        figureParams;
    end
    properties (Dependent)
        currTime;
    end
    properties (Hidden)
        CDlength;
        isExpClosed=0;
        isTraining=0;
    end
    methods
        % Methods
        function obj=initialize(obj,varargin)
            % If an argument is passed, it must be a structure matching
            % linearMap template (i.e. basically, an output of another
            % session). Some parameters (e.g. sampling frequency of
            % amplifier and buffer window length) cannot be changed
            % directly from here, so please do make sure that they're
            % matching with the current settings of relevant Simulink
            % model.
            
            % Normalization buffer length
            obj.bufferData.bufferLength=30; % Data will be normalized in mean and SD over this amount of seconds
            
            % Default parameters
            if isempty(varargin{1})
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
            % Variables on base workspace will be used to trigger closing
            % of experiment
            assignin('base','isExpClosing',0);
            
            % Sets name of Simulink model to be used for acquisition
            obj.modelName='SimpleAcquisition_16ch_2014a_RT';
            
            % Prepares Simulink model (i.e. starts recording, basically)
            obj.prepareSimulinkModel;
            
            % Generates array of time triggered events
            obj.timeTriggeredEvents{1}=timeTriggeredEvent('expCallback',0);
            obj.timeTriggeredEvents{2}=timeTriggeredEvent('toggleTraining',0);
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
            
            % Clear variables from base workspace
            evalin('base','clear listener*');
            evalin('base','clear toggleTraining');
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
            obj.linearMap.mat=MI_ET.updateWeights(obj.linearMap.mat,BP,currEst,obj.inTarget,1e-3);
        end
        
        function obj=toggleTraining(obj)
            if evalin('base','exist(''toggleTraining'',''var'')')&&evalin('base','toggleTraining')
                obj.isTraining=~obj.isTraining;
                assignin('base','toggleTraining',0);
                figure(obj.figureParams.handle)
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
        
        %% Dependent properties
        function cTime=get.currTime(obj)
            if obj.isExpClosed
                cTime=obj.rawData.Time(end);
            else
                cTime=get_param(obj.modelName,'SimulationTime');
            end
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
                MI_exp.closeExp;
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
            MI_exp.closeExp;
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