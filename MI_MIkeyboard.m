classdef MI_MIkeyboard
    properties
        fs=512; % WARNING: DO NOT change this. Amplifier acquisition rate can only be changed from its panel in Simulink, it is here for reference
        rawData;
        condition;
        modelName;
        timeTriggeredEvents;
        outputLog;
        timingParams;
        keyboardUDPchannels;
        mat18UDPchannels;
        lastKeyID;
        selStartTime;
        selStartID='';
        selCounter=0;
        selDuration;
        selDefaultDuration=2;
        selSound;
        keyPressTimeBuffer=[];
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
        UDPtimeOut=0.2;
    end
    methods
        %% Constructor
        function obj=MI_MIkeyboard
            % Some parameters (e.g. sampling frequency of amplifier and
            % buffer window length) cannot be changed directly from here,
            % so please do make sure that they're matching with the current
            % settings of relevant Simulink model.
            
            % Set length of initial countdown, in seconds
            obj.CDlength=15;
            
            % Define expected inputs and output ports for udp communication
            obj.keyboardUDPchannels.outPort=9000;
            obj.keyboardUDPchannels.inPort=8051;
            obj.mat18UDPchannels.outPort=7561;
            obj.mat18UDPchannels.inPort=7560;
            
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
            obj.outputLog.keyUDPlogIn=[];
            obj.outputLog.keyUDPlogOut=[];
            % Following key selection
            obj.outputLog.selStart=[];
            obj.outputLog.selEnd=[];
            obj.outputLog.selID=[];
            obj.outputLog.selDurationList=[];
            
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
            obj.mat18UDPchannels.udps.step(uint8('close'));
            save(fileName,'obj');
            
            % Release receiver UDP port
            obj.keyboardUDPchannels.udpr.release;
            obj.mat18UDPchannels.udpr.release;

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
            obj.keyboardUDPchannels.lastInput=char(obj.keyboardUDPchannels.udpr.step)';
            inputTime=obj.currTime;
            
            % Perform the following if gaze is moving into a key
            if ~strcmp(obj.selStartID,obj.currKeyID)||(obj.outputLog.isInTarget(end)==0&&obj.inTarget==1)
                % Recover most recent raw data buffer from workspace and
                % send it to Matlab2018 over UDP for analysis
                dataPacket=evalin('base','bufferData');
                outString=uint8(mat2str(dataPacket,4));
                obj.mat18UDPchannels.udps.step(outString);
                waitStartTime=obj.currTime;
                
                % Wait briefly for net response from UDP. If there is no
                % response, move on with default duration
                obj.selDuration=obj.selDefaultDuration;
                while obj.currTime-waitStartTime<obj.UDPtimeOut
                    %                     fprintf('%0.2f, %0.2f\n',obj.currTime-waitStartTime,obj.UDPtimeOut)
                    wait(obj,0.01);
                    inString=obj.mat18UDPchannels.udpr.step;
                    if ~isempty(inString)
                        obj.selDuration=eval(char(inString));
                        break;
                    end
                end
                fprintf('%0.1f\n',obj.selDuration)
                
                obj.selStartID=obj.currKeyID;
                obj.selStartTime=inputTime;
                obj.selCounter=0;
            end
            
            % If cursor is within target, fills up selection counter
            if obj.inTarget
                obj.selCounter=(obj.currTime-obj.selStartTime)/obj.selDuration;
            end
                                                                                                                  
            % If current key has already been selected, reset counter to zero
            if obj.selCounter>=1
                % Play sound
                sound(obj.selSound.y,obj.selSound.f);
                
                % Communicate whether threshold has been reached
                outString=uint8(sprintf('%d',1));
                
                % Reset selection counter and update relevant logs
                obj.selCounter=0;
                obj.outputLog.selStart=cat(1,obj.outputLog.selStart,obj.selStartTime);
                obj.outputLog.selEnd=cat(1,obj.outputLog.selEnd,inputTime);
                obj.outputLog.selID{end+1}=obj.currKeyID;
                obj.outputLog.selDurationList=cat(1,obj.outputLog.selDurationList,obj.selDuration);
                obj.selStartTime=obj.currTime;
                
                % Update buffer of key presses
                obj.keyPressTimeBuffer=[obj.keyPressTimeBuffer,obj.currTime];
            else
                outString=uint8(sprintf('%0.2f',min(obj.selCounter,0.989))); % Threshold output to prevent sending 1 as result of a rounding error
            end
            obj.keyboardUDPchannels.udps.step(outString);
            
            % Add relevant info to log
            obj.outputLog.time=cat(1,obj.outputLog.time,obj.currTime);
            obj.outputLog.isTraining=cat(1,obj.outputLog.isTraining,obj.isTraining);
            obj.outputLog.isInTarget=cat(1,obj.outputLog.isInTarget,obj.inTarget);
            obj.outputLog.selCounter=cat(1,obj.outputLog.selCounter,obj.selCounter);
            obj.outputLog.keyUDPlogIn{end+1}=obj.keyboardUDPchannels.lastInput;
            obj.outputLog.keyUDPlogOut{end+1}=outString;
            obj.outputLog.cursorPos=cat(1,obj.outputLog.cursorPos,obj.ETcursorPos);
            
            % Set next evaluation time for this function
            obj.timeTriggeredEvents{1}.triggersLog=[obj.timeTriggeredEvents{1}.triggersLog,obj.currTime];
            obj.timeTriggeredEvents{1}.nextTrigger=obj.currTime;
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
            obj.keyboardUDPchannels.udpr=dsp.UDPReceiver('LocalIPPort',obj.keyboardUDPchannels.inPort,'ReceiveBufferSize',1);
            obj.keyboardUDPchannels.udps=dsp.UDPSender('RemoteIPPort',obj.keyboardUDPchannels.outPort);
            obj.mat18UDPchannels.udpr=dsp.UDPReceiver('LocalIPPort',obj.mat18UDPchannels.inPort,'ReceiveBufferSize',1);
            obj.mat18UDPchannels.udps=dsp.UDPSender('RemoteIPPort',obj.mat18UDPchannels.outPort);
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
 
        function obj=attachPhase(obj,otherSession)
            timeStep=median(diff(obj.rawData.time));
            otherSession.rawData.time=linspace(obj.rawData.time(end)+timeStep,obj.rawData.time(end)+otherSession.rawData.time(end),length(otherSession.rawData.time));
            joiningFields={'cursorPos','selCounter','isTraining','isInTarget','keyUDPlogOut','keyUDPlogIn'};
            updatingFields={'time','selStart','selEnd'};
            for currUpField=1:length(updatingFields)
                obj.outputLog.(updatingFields{currUpField})=cat(1,obj.outputLog.(updatingFields{currUpField}),otherSession.outputLog.(updatingFields{currUpField})+obj.rawData.time(end));
            end
            for currJoinField=1:length(joiningFields)
                try
                obj.outputLog.(joiningFields{currJoinField})=cat(find(size(obj.outputLog.(joiningFields{currJoinField}))==max(size(obj.outputLog.(joiningFields{currJoinField})))),obj.outputLog.(joiningFields{currJoinField}),otherSession.outputLog.(joiningFields{currJoinField}));
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
            if isfield(obj.keyboardUDPchannels,'lastInput')&&~isempty(obj.keyboardUDPchannels.lastInput)&&~sum(isnan(obj.keyboardUDPchannels.lastInput))
                if numel(strfind(obj.keyboardUDPchannels.lastInput,'True'))
                    inT=1;
                end
            end
        end
        
        function keyID=get.currKeyID(obj)
            if obj.inTarget&&isfield(obj.keyboardUDPchannels,'lastInput')&&~isempty(obj.keyboardUDPchannels.lastInput)&&~sum(isnan(obj.keyboardUDPchannels.lastInput))
                semiColumnPos=strfind(obj.keyboardUDPchannels.lastInput,';');
                keyID=obj.keyboardUDPchannels.lastInput(semiColumnPos(4)+1:semiColumnPos(5));
            else
                keyID=0;
            end
%             fprintf('%d\n',keyID)
        end
        
        function pos=get.ETcursorPos(obj)
            if isfield(obj.keyboardUDPchannels,'lastInput')&&~isempty(obj.keyboardUDPchannels.lastInput)&&~sum(isnan(obj.keyboardUDPchannels.lastInput))
                lastEntry=obj.keyboardUDPchannels.lastInput;
                semiColumnPos=strfind(lastEntry,';');
                lastEntry(semiColumnPos)=',';
                pos=str2num(lastEntry(2:semiColumnPos(2)-1)); %#ok<ST2NM>
            else
                pos=nan(1,2);
            end
        end
    end
    methods (Static)
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
blockName=sprintf('%s/bufferData/rawData_buffer',modelName);
assignin('base','listenerMI',add_exec_event_listener(blockName,'PostOutputs',@acquireBufferData));
end

function acquireBufferData(block,~)
assignin('base','bufferData',block.OutputPort(1).Data);
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
    MI_MIkeyboard.closeExp;
end
if strcmp(eventdata.Key,'p')
    keyboard;
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
MI_MIkeyboard.closeExp;
end