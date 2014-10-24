%**************************************************************************
% Class definition for a utility to allow the user to display a struct or
% class hierarchy.
%
% Usage:
%         % instantiate the utility; figure appears and is ready to have
%         % have the hieratchy drawn on it
%         >> util = InstanceHierarchyUtility; 
%
%         % pass an any struct or object to the util to display its
%         % hierarchy on the plot
%         >> util.show(x);
%
% Notes:
%        Developed and tested for MATLAB R2014a.
%
%**************************************************************************
% Code is written by author based upon noted references (if given), using 
% The MathWorks MATLAB function signature (if applicable) for convenience 
% only.
%
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
classdef InstanceHierarchyUtility < handle
    
    % private data members
    properties (Access=private)       
        % Handle to the figure displaying the hierarchy.
        figHandle;
        % Structure or class instance to display the hierarchy for.
        instance;
    end
    
    % public data members
    properties (Access=public)
        % intentionally left blank
    end
    
    % events for anyone listening
    events
    end
    
    methods (Access=public)
       
        % Utility class constructor.
        %
        % @param   inst   Structure or class instance to draw the hierarchy
        %                 for.
        % @param   levels The number of levels deep in the hierarchy in the
        %                 to display.
        %
        % @return  The plot markup utility
        %
        function [util] = InstanceHierarchyUtility()
            
            util = init(util);
            
            util = util.show;
        end 
        
        % Removes the last markup to the plot.
        %
        % @return  The plot markup utility
        %
        function [util] = removeLastMarkup(util)
            
            if ~isempty(util.markupHandles)
                delete(util.markupHandles(end));
                util.markupHandles = util.markupHandles(1:end-1);
            end
        end 
        
        % Set the markup colour
        %
        % @param   colour   The new markup colour.
        %
        % @return  The plot markup utility
        %
        function [util] = setMarkupColour(util,colour)            
            util = validateColour(util,colour);
        end         
        
        % Set the markup width
        %
        % @param   width   The new markup line width.
        %
        % @return  The plot markup utility
        %
        function [util] = setMarkupLineWidth(util,width)            
            util.markupLineWidth = max(0.5,min(width,50));
        end          
        
        % Closes the existing figure.
        %
        % @return  The plot markup utility
        %
        function [util] = closeFigure(util)
            util = close(util);
        end    
        
        % Creates a new figure.
        %
        % @return  The plot markup utility
        %
        function [util] = newFigure(util)
            newFigPosition = [];
            if ~isempty(util.figHandle)
                newFigPosition = get(util.figHandle,'Position');
            end
            util = close(util);
            util = show(util);
            
            if ~isempty(newFigPosition)
                set(util.figHandle,'Position',newFigPosition);
            end            
        end          
        
        
        % Plot markup utility class destructor.
        %
        % @param   util   The plot markup utility instance to 
        %                         delete.
        %
        % @return  The plot markup utility instance.
        %
        function [util] = delete(util)
            util = close(util);
            util = init(util);
        end
        
        % Displays context for the object.
        %
        % @param   util   The plot markup utility to describe.
        %
        function display(util)
            fprintf('\n%s = \n\t',inputname(1));
            fprintf('isButtonDown      %d\n\t',util.isButtonDown);
            if isempty(util.figHandle)
                fprintf('figHandle     []\n\t');
            else
                fprintf('figHandle     %f\n\t', util.figHandle);
            end
            if isempty(util.plotAxisHandle)
                fprintf('plotAxisHandle    []\n\t');
            else
                fprintf('plotAxisHandle    %f\n\t', util.plotAxisHandle);
            end            
            if isempty(util.markupHandles)
                fprintf('markupHandles     []\n\t');
            else
                fprintf('markupHandles     [%dx%d]\n\t', ...
                    size(util.markupHandles,1),size(util.markupHandles,2));
            end   
            fprintf('markupColour      [%d %d %d]\n\t',util.markupColour);
            fprintf('markupLineWidth   %f\n',util.markupLineWidth);
            fprintf('\n');
        end
        
        % Callback for the mouse down event.
        function mouseDownCallback(util,~,~)

            mouseDownPoint = ...
                get(util.plotAxisHandle, 'CurrentPoint');
            mouseDownPoint = mouseDownPoint(1,1:2);

            % only record event if the button is down within the image
            xlim = get(util.plotAxisHandle,'XLim');
            ylim = get(util.plotAxisHandle,'YLim');

            if mouseDownPoint(1)>=xlim(1) &&  ...
               mouseDownPoint(1)<=xlim(2) &&  ...
               mouseDownPoint(2)>=ylim(1) &&  ...
               mouseDownPoint(2)<=ylim(2)
                util.isButtonDown     = true;
                util.firstMarkupPoint = mouseDownPoint';
                hold on;
            else
                util.isButtonDown     = false;
                util.markupHandles    = [];
                util.firstMarkupPoint = [];
            end

        end

        % Callback for the mouse motion event.
        function mouseMotionCallback(util,~,~)

            if util.isButtonDown

                mouseMotionPoint = ...
                    get(util.plotAxisHandle, 'CurrentPoint');
                mouseMotionPoint = mouseMotionPoint(1,1:2);
                
                % only record event if mouse pointer is within axes
                xlim = get(util.plotAxisHandle,'XLim');
                ylim = get(util.plotAxisHandle,'YLim');

                if mouseMotionPoint(1)>=xlim(1) &&  ...
                   mouseMotionPoint(1)<=xlim(2) &&  ...
                   mouseMotionPoint(2)>=ylim(1) &&  ...
                   mouseMotionPoint(2)<=ylim(2)
               
                    if isempty(util.firstMarkupPoint)
                        hndle = util.markupHandles(end);
                        xData = [get(hndle,'XData') mouseMotionPoint(1)];
                        yData = [get(hndle,'YData') mouseMotionPoint(2)];
                        set(hndle,'XData',xData,'YData',yData);
                    else
                        origin = util.firstMarkupPoint;
                        h      = plot([origin(1);mouseMotionPoint(1)],[origin(2);mouseMotionPoint(2)],...
                            '-','Color',util.markupColour,'LineWidth',util.markupLineWidth);
                        
                        util.markupHandles    = [util.markupHandles; h];
                        util.firstMarkupPoint = [];
                    end
                end  
            end
        end

        % Callback for the mouse button up event.
        function mouseUpCallback(util,~,~)

            if util.isButtonDown
                
                util.isButtonDown = false;
                
                % Broadcast notice of event
                %notify(util,'PAUMouseButtonUp',IRUMouseButtonUpEvent(util.imgRoi)); 
                
            end  
        end        
    end
    
    methods (Access=private)
        
        % Displays the image in a figure.
        %
        % @param   util   The plot markup utility to show the
        %                           figure for.
        %
        function [util] = show(util)
            
            % create the figure
            util.figHandle  = figure;
            util.plotAxisHandle = axes;

            % set the callback functions for the figure
            set(util.figHandle, ...
                'WindowButtonDownFcn',   @util.mouseDownCallback, ...
                'WindowButtonUpFcn',     @util.mouseUpCallback,   ...
                'WindowButtonMotionFcn', @util.mouseMotionCallback);    
        end
        
        % Initializes the data members of the plot markup utility.
        %
        % @param   util   The plot markup instance to initialize.
        %
        % @return  The plot markup utility instance.
        %       
        function [util] = init(util)
            util.isButtonDown     = false;
            util.figHandle    = [];
            util.plotAxisHandle   = [];
            util.markupHandles    = [];
            util.firstMarkupPoint = [];
            util.markupColour     = [1 0 0];
            util.markupLineWidth      = 0.5;
        end
        
        % Closes the currently open figure.
        %
        % @param   util   The plot markup instance to initialize.
        %
        % @return  The plot markup utility instance.
        %       
        function [util] = close(util)
            if ~isempty(util.figHandle)                
                
                for i=1:length(util.markupHandles)
                    delete(util.markupHandles(i));
                end
                
                util.markupHandles = [];
                
                try
                    close(util.figHandle);
                    util.figHandle  = [];
                    util.plotAxisHandle = [];
                catch
                    % intentionally left blank
                end
            end
        end        
        
        % Validates the markup colour.
        %
        % @param   colour   The colour to validate.
        %
        % @return  The plot markup utility instance.
        function [util] = validateColour(util,colour)
            
            if ischar(colour)
                colour = lower(colour);
                switch colour
                    case {'y','yellow'}
                        util.markupColour = [1 1 0];
                    case {'m','magenta'}
                        util.markupColour = [1 0 1];
                    case {'c','cyan'}
                        util.markupColour = [0 1 1];
                    case {'r','red'}
                        util.markupColour = [1 0 0];
                    case {'g','green'}
                        util.markupColour = [0 1 0];
                    case {'b','blue'}
                        util.markupColour = [0 0 1];
                    case {'w','white'}
                        util.markupColour = [1 1 1];
                    case {'k','black'}
                        util.markupColour = [0 0 0];
                    otherwise
                        printf('InstanceHierarchyUtility - colour string is invalid');
                end
            else
                if length(colour)~=3
                    fprintf('InstanceHierarchyUtility - colour vector must have 3 elements');
                elseif ~all(colour<=1 & colour>=0.0)
                    fprintf('InstanceHierarchyUtility - colour vector elements must be between 0 and 1');
                else 
                    util.markupColour = double(colour);
                end
            end
        end
    end
end

