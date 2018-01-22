classdef MI_ET_sub < MI_exp
    properties
        cursorPos;
        targetPos;
        ETdata;
        recLength
        nTargets;
        currTrial;
        actualTarget;
        inTarget;
        selCounter=0;
        selThreshold=.9;
        selDecay=.1;
    end
    properties (Hidden)
        gazePosBuffer; % Buffer of gaze positions within each target
        tform; % Calibration transformation from pixel coordinates of gaze to image coordinates        
    end
    methods
        %% Constructor
        function obj=MI_ET_sub(varargin)
            % If an argument is passed, it must be a structure matching
            % linearMap template (i.e. basically, an output of another
            % session). Some parameters (e.g. sampling frequency of
            % amplifier and buffer window length) cannot be changed
            % directly from here, so please do make sure that they're
            % matching with the current settings of relevant Simulink
            % model.
            
            % Set desired length of recording. 
            obj.recLength=600;
                        
            % Set length of initial countdown, in seconds
            obj.CDlength=15;
            
            % Set colors for different objects
            obj.figureParams.bg=[.05,.05,.05];
            obj.figureParams.targetColor=[0,.4,0];
            obj.figureParams.cursorColor=[.4,0,.1];
            
            % Set shape and pos for cursor
            obj.figureParams.cursorShape.X=[-.05,.05,.05,-.05];
            obj.figureParams.cursorShape.Y=[-.05,-.05,.05,.05];
            obj.cursorPos=[0,0];
            
            % Set possible target centers and shape
            targetX=(1:3)'*ones(1,3);
            targetY=ones(3,1)*(1:3);
            targetX=(targetX-2)*.9;
            targetY=(targetY-2)*.9;
            obj.nTargets=numel(targetX);
            obj.targetPos=cell(obj.nTargets,1);
            for currT=1:obj.nTargets
                obj.targetPos{currT}=[targetX(currT),targetY(currT)];
            end
            obj.figureParams.targetRadius=.1;
            obj.figureParams.targetShape.X=cos(linspace(0,2*pi,40))*obj.figureParams.targetRadius;
            obj.figureParams.targetShape.Y=sin(linspace(0,2*pi,40))*obj.figureParams.targetRadius;
            
            % Perform initialization
            obj=initialize(obj,varargin);
                        
            % Initialize a few more things
            obj.outputLog.cursorPos=[];
            obj.outputLog.targetsReached.time=[];
            obj.outputLog.actualTarget=[];
            obj.outputLog.isInTarget=[];
            obj.outputLog.selCounter=[];
            obj.outputLog.currEst=[];
            
            % Ask user whether to start experiment right away
            clc;
            if ~strcmpi(input('Start experiment now? [Y/n]\n','s'),'n')
                obj=start(obj);
            end
        end
        
        % Other methods
        function obj=start(obj)            
            % Start eye-tracking
            MI_ET_sub.startEyeTracking;
            
            % Ask user whether to calibrate eye-tracking. Cursor should not
            % be visible during calibration
            persistent calibrationTransform
            if ~isempty(calibrationTransform)
                userChoice=input('Perform gaze tracking calibration? [y/N]\n','s');
            else
                userChoice='y';
            end
                        
            % Opens black figure as background
            obj=createExpFigure(obj);
            
            % Perform calibration, if needed
            if strcmpi(userChoice,'y');
                obj=obj.performETcalibration;
                calibrationTransform=obj.tform;
            end
            
            % Randomly select a target
            obj.actualTarget=ceil(rand*obj.nTargets);
            
            % Shows a countdown
            obj.startCountdown(obj.CDlength);
            
            % Perform bulk of experiment
            obj=manageExperiment(obj);
            
            % Closes exp window and saves data
            obj.closeExp;
        end
                
        function obj=createExpFigure(obj)
            % Set figure properties
            obj.figureParams.handle=gcf;
            set(obj.figureParams.handle,'Tag',mfilename,...
                'Toolbar','none',...
                'MenuBar','none',...
                'Units','normalized',...
                'Resize','off',...
                'NumberTitle','off',...
                'Name','',...
                'Color',obj.figureParams.bg,...
                'RendererMode','Manual',...
                'Renderer','OpenGL',...
                'WindowKeyPressFcn',@MI_exp.KeyPressed,...
                'CloseRequestFcn',@MI_exp.OnClosing,...
                'WindowButtonMotionFcn',@MI_exp.onMouseMove);
            
            % Create cursor outside of visible field
            obj.figureParams.cursor=patch(obj.figureParams.cursorShape.X+100,obj.figureParams.cursorShape.Y+100,obj.figureParams.cursorColor);
            
            % Set and remove figure axis
            ylim([-1,1]);
            xlim([-1,1]);
            set(gcf,'units','normalized','position',[0,0,1,1]);
            axis square
            axis('off')
            
            % Remove box around figure
            %             undecorateFig;
        end
        
        function obj=expCallback(obj)
            global udpr
            % Evaluate cursor position from gaze
            coords=str2num(char(udpr.step)'); %#ok<ST2NM>
            if ~isempty(coords)
                obj.cursorPos=obj.tform.transformPointsForward(coords');
            end
            
            % Filter gaze coordinates
            persistent previousGazePos
            if isempty(previousGazePos)||sum(isnan(previousGazePos))
                previousGazePos=obj.cursorPos;
            end
            gazeSpeed=sqrt(sum((previousGazePos-obj.cursorPos).^2));
            obj.cursorPos=previousGazePos*(1-sqrt(gazeSpeed))+obj.cursorPos*sqrt(gazeSpeed);
            previousGazePos=obj.cursorPos;
            set(obj.figureParams.cursor,'XData',obj.figureParams.cursorShape.X+obj.cursorPos(1),'YData',obj.figureParams.cursorShape.Y+obj.cursorPos(2));
            
            % Check if cursor is within target
            if sqrt(sum((obj.cursorPos-obj.targetPos{obj.actualTarget}).^2))<obj.figureParams.targetRadius
                obj.inTarget=1;
            else
                obj.inTarget=0;
            end
            
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
                [~,obj.linearMap.lapFltrWeights]=MI_ET_sub.applyLapFilter(dataWindow);
            end
            
            % Recover bandpower data from data buffer
            BP=MI_ET_sub.preprocData(dataWindow,obj.linearMap);
            
            % If training is ongoing, update linear map
            if obj.isTraining
                obj=obj.computeLinearMap(BP);
            end
            
            % If cursor is within target, fills up selection counter
            if isfield(obj.linearMap,'mat')
                feats=cat(2,1,BP);
                currEst=1./(1+exp(-(feats*obj.linearMap.mat)));
            else
                currEst=0;
            end
            if obj.inTarget
                % Use EEG data when available, otherwise use gaze alone
                if isfield(obj.linearMap,'mat')
                    obj.selCounter=max(0,obj.selCounter+currEst-obj.selThreshold);
                else
                    obj.selCounter=obj.selCounter+.01;
                end
                if obj.selCounter>=1
                    obj.actualTarget=ceil(rand*obj.nTargets);
                    obj.selCounter=0;
                end
%                 
%                 % If obj is within a target, update gaze positions buffer
%                 obj.gazePosBuffer{obj.actualTarget}=cat(2,obj.gazePosBuffer{obj.actualTarget},obj.cursorPos');
%                 obj.gazePosBuffer{obj.actualTarget}(:,1)=[];
%                 obj=obj.calibrateEyeTracker;
            else
                obj.selCounter=obj.selCounter*(1-obj.selDecay);
            end
            
            % Update target pos and color
            xTarget=obj.figureParams.targetShape.X+obj.targetPos{obj.actualTarget}(1);
            yTarget=obj.figureParams.targetShape.Y+obj.targetPos{obj.actualTarget}(2);
            targetColor=obj.figureParams.targetColor*(1-obj.selCounter);
            set(obj.figureParams.target,'XData',xTarget,'YData',yTarget,'FaceColor',targetColor);
            drawnow;
            
            % Add relevant info to log
            obj.outputLog.cursorPos=cat(1,obj.outputLog.cursorPos,obj.cursorPos);
            obj.outputLog.actualTarget=cat(1,obj.outputLog.actualTarget,obj.actualTarget);
            obj.outputLog.feats=cat(1,obj.outputLog.feats,BP);
            obj.outputLog.isTraining=cat(1,obj.outputLog.isTraining,obj.isTraining);
            obj.outputLog.time=cat(1,obj.outputLog.time,obj.currTime);
            obj.outputLog.isInTarget=cat(1,obj.outputLog.isInTarget,obj.inTarget);
            obj.outputLog.selCounter=cat(1,obj.outputLog.selCounter,obj.selCounter);
            obj.outputLog.currEst=cat(1,obj.outputLog.currEst,currEst);
            
            % Set next evaluation time for this function
            obj.timeTriggeredEvents{1}.triggersLog=[obj.timeTriggeredEvents{1}.triggersLog,obj.currTime];
            obj.timeTriggeredEvents{1}.nextTrigger=obj.currTime+.05;
        end
                                
        function obj=performETcalibration(obj)
            % Generate one target in each position for about 200 steps
            % (each step is 10 ms).
            global udpr
            obj.figureParams.target=patch(obj.figureParams.targetShape.X,obj.figureParams.targetShape.Y,obj.figureParams.targetColor);
            obj.gazePosBuffer=cell(obj.nTargets,1);
            pause(2);
            for currTarget=1:obj.nTargets
                set(obj.figureParams.target,'XData',obj.targetPos{currTarget}(1)+obj.figureParams.targetShape.X,'YData',obj.targetPos{currTarget}(2)+obj.figureParams.targetShape.Y);
                drawnow;
                obj.gazePosBuffer{currTarget}=[];
                for currRep=1:200
                    coords=str2num(char(udpr.step)'); %#ok<ST2NM>
                    obj.gazePosBuffer{currTarget}=cat(2,obj.gazePosBuffer{currTarget},coords);
                    pause(0.01);
                end
            end
            
            % Move target off-screen
            set(obj.figureParams.target,'XData',100+obj.figureParams.targetShape.X,'YData',obj.targetPos{currTarget}(2)+obj.figureParams.targetShape.Y);
            
            % Perform actual calibration
            obj=obj.calibrateEyeTracker;
        end
       
        function obj=calibrateEyeTracker(obj)
            % Use coordinates of gaze during each target fixation to
            % calibrate transformation
            fixedPoints=cell2mat(obj.targetPos);
            fixedPointsMat=[];
            for currLbl=1:9
                fixedPointsMat=cat(1,fixedPointsMat,repmat(fixedPoints(currLbl,:),length(obj.gazePosBuffer{currLbl}),1));
            end
            movingPointsMat=cell2mat(obj.gazePosBuffer')';
            obj.tform=fitgeotrans(movingPointsMat,fixedPointsMat,'Projective');
        end
                
        function startCountdown(obj,nSecs)
            % countdown to experiment start
            figure(obj.figureParams.handle)
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
    end
    methods (Static)                
        function startEyeTracking
            % Losing udpr is bad (cannot close port from Matlab anymore).
            % Declare it as global so that it can be recovered outside of
            % this class, as well
            global udpr
            % Local port numeber is given by GazeTrackEyeXGazeStream
            % Buffer size is pretty much arbitrarily chosen: it should be
            % so that it contains little more than one entry
            if isempty(udpr)
                udpr=dsp.UDPReceiver('LocalIPPort',11000,'ReceiveBufferSize',20);
            end
            
            % Launch GazeTrackEyeXGazeStream in async mode, if not already
            % running (prompts user, I have no idea how to check if
            % external code is running)
            clc
            startET=input('WARNING: select Yes only if GazeStream is not running already.\nDo you want to start eye tracking? y/N: ','s');
            if isempty(startET)
                startET='n';
            end
            if strcmpi(startET,'y')
                !C:\Code\Sources\GazeTrackEyeXGazeStream\GazeTrackEyeXGazeStream.exe &
            end
        end
    end
end