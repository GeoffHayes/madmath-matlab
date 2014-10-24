%**************************************************************************
% Class definition for a utility to allow the user to mark up a plot.
%
% Usage:
%         % instantiate the utility; figure appears and is ready to have
%         % have plots/lines/shapes etc. drawn on it
%         >> plotMarkupUtil = PlotMarkupUtility; 
%
%         % run any command to draw something on the figure
%         >> x=-5:01:5;y=x.^2;
%         >> plot(x,y);
%
%         % move the cursor within the figure axes and press the mouse
%         % button
%
%         % while keeping the button down, move the cursor which leaves a
%         % trail of points corresponding to the movement of the mouse
%
%         % release mouse button when the markup is complete
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
classdef PlotMarkupUtility < handle
    
    % private data members
    properties (Access=private)       
        % Indicator as to whether the user has pressed the mouse button
        % within the image figure.
        isButtonDown;
        % Handle to the figure displaying the plot.
        plotFigHandle;
        % Handle to the axis displaying the plot.
        plotAxisHandle;        
        % Set of handles corresponding to points that have been used to
        % markup the plot (one handles for every set of points).
        markupHandles;
        % First point recorded on a mouse down event
        firstMarkupPoint;
        % The markup colour.
        markupColour;
        % The markup line width.
        markupLineWidth;
    end
    
    % public data members
    properties (Access=public)
        % intentionally left blank
    end
    
    % events for anyone listening
    events
      PAUMouseButtonUp % fired on a mouse button up/release
    end
    
    methods (Access=public)
       
        % Plot markup utility class constructor.
        %
        % @param   colour   Markup colour.
        %
        % @return  The plot markup utility
        %
        function [plotMarkupUtil] = PlotMarkupUtility(colour)
            
            plotMarkupUtil = init(plotMarkupUtil);
            
            if exist('colour','var')
                plotMarkupUtil = validateColour(plotMarkupUtil,colour);
            end
            
            plotMarkupUtil = plotMarkupUtil.show;
        end 
        
        % Removes the last markup to the plot.
        %
        % @return  The plot markup utility
        %
        function [plotMarkupUtil] = removeLastMarkup(plotMarkupUtil)
            
            if ~isempty(plotMarkupUtil.markupHandles)
                delete(plotMarkupUtil.markupHandles(end));
                plotMarkupUtil.markupHandles = plotMarkupUtil.markupHandles(1:end-1);
            end
        end 
        
        % Set the markup colour
        %
        % @param   colour   The new markup colour.
        %
        % @return  The plot markup utility
        %
        function [plotMarkupUtil] = setMarkupColour(plotMarkupUtil,colour)            
            plotMarkupUtil = validateColour(plotMarkupUtil,colour);
        end         
        
        % Set the markup width
        %
        % @param   width   The new markup line width.
        %
        % @return  The plot markup utility
        %
        function [plotMarkupUtil] = setMarkupLineWidth(plotMarkupUtil,width)            
            plotMarkupUtil.markupLineWidth = max(0.5,min(width,50));
        end          
        
        % Closes the existing figure.
        %
        % @return  The plot markup utility
        %
        function [plotMarkupUtil] = closeFigure(plotMarkupUtil)
            plotMarkupUtil = close(plotMarkupUtil);
        end    
        
        % Creates a new figure.
        %
        % @return  The plot markup utility
        %
        function [plotMarkupUtil] = newFigure(plotMarkupUtil)
            newFigPosition = [];
            if ~isempty(plotMarkupUtil.plotFigHandle)
                newFigPosition = get(plotMarkupUtil.plotFigHandle,'Position');
            end
            plotMarkupUtil = close(plotMarkupUtil);
            plotMarkupUtil = show(plotMarkupUtil);
            
            if ~isempty(newFigPosition)
                set(plotMarkupUtil.plotFigHandle,'Position',newFigPosition);
            end            
        end          
        
        
        % Plot markup utility class destructor.
        %
        % @param   plotMarkupUtil   The plot markup utility instance to 
        %                         delete.
        %
        % @return  The plot markup utility instance.
        %
        function [plotMarkupUtil] = delete(plotMarkupUtil)
            plotMarkupUtil = close(plotMarkupUtil);
            plotMarkupUtil = init(plotMarkupUtil);
        end
        
        % Displays context for the object.
        %
        % @param   plotMarkupUtil   The plot markup utility to describe.
        %
        function display(plotMarkupUtil)
            fprintf('\n%s = \n\t',inputname(1));
            fprintf('isButtonDown      %d\n\t',plotMarkupUtil.isButtonDown);
            if isempty(plotMarkupUtil.plotFigHandle)
                fprintf('plotFigHandle     []\n\t');
            else
                fprintf('plotFigHandle     %f\n\t', plotMarkupUtil.plotFigHandle);
            end
            if isempty(plotMarkupUtil.plotAxisHandle)
                fprintf('plotAxisHandle    []\n\t');
            else
                fprintf('plotAxisHandle    %f\n\t', plotMarkupUtil.plotAxisHandle);
            end            
            if isempty(plotMarkupUtil.markupHandles)
                fprintf('markupHandles     []\n\t');
            else
                fprintf('markupHandles     [%dx%d]\n\t', ...
                    size(plotMarkupUtil.markupHandles,1),size(plotMarkupUtil.markupHandles,2));
            end   
            if isempty(plotMarkupUtil.firstMarkupPoint)
                fprintf('firstMarkupPoint  []\n\t');
            else
                fprintf('firstMarkupPoint  [%f,%f]\n\t', ...
                    plotMarkupUtil.firstMarkupPoint(1),plotMarkupUtil.firstMarkupPoint(2));
            end
            fprintf('markupColour      [%d %d %d]\n\t',plotMarkupUtil.markupColour);
            fprintf('markupLineWidth   %f\n',plotMarkupUtil.markupLineWidth);
            fprintf('\n');
        end
        
        % Callback for the mouse down event.
        function mouseDownCallback(plotMarkupUtil,~,~)

            mouseDownPoint = ...
                get(plotMarkupUtil.plotAxisHandle, 'CurrentPoint');
            mouseDownPoint = mouseDownPoint(1,1:2);

            % only record event if the button is down within the image
            xlim = get(plotMarkupUtil.plotAxisHandle,'XLim');
            ylim = get(plotMarkupUtil.plotAxisHandle,'YLim');

            if mouseDownPoint(1)>=xlim(1) &&  ...
               mouseDownPoint(1)<=xlim(2) &&  ...
               mouseDownPoint(2)>=ylim(1) &&  ...
               mouseDownPoint(2)<=ylim(2)
                plotMarkupUtil.isButtonDown     = true;
                plotMarkupUtil.firstMarkupPoint = mouseDownPoint';
                hold on;
            else
                plotMarkupUtil.isButtonDown     = false;
                plotMarkupUtil.markupHandles    = [];
                plotMarkupUtil.firstMarkupPoint = [];
            end

        end

        % Callback for the mouse motion event.
        function mouseMotionCallback(plotMarkupUtil,~,~)

            if plotMarkupUtil.isButtonDown
                mouseMotionPoint = ...
                    get(plotMarkupUtil.plotAxisHandle, 'CurrentPoint');
                mouseMotionPoint = mouseMotionPoint(1,1:2);
                
                % only record event if mouse pointer is within axes
                xlim = get(plotMarkupUtil.plotAxisHandle,'XLim');
                ylim = get(plotMarkupUtil.plotAxisHandle,'YLim');

                if mouseMotionPoint(1)>=xlim(1) &&  ...
                   mouseMotionPoint(1)<=xlim(2) &&  ...
                   mouseMotionPoint(2)>=ylim(1) &&  ...
                   mouseMotionPoint(2)<=ylim(2)
               
                    if isempty(plotMarkupUtil.firstMarkupPoint)
                        hndle = plotMarkupUtil.markupHandles(end);
                        xData = [get(hndle,'XData') mouseMotionPoint(1)];
                        yData = [get(hndle,'YData') mouseMotionPoint(2)];
                        set(hndle,'XData',xData,'YData',yData);
                    else
                        origin = plotMarkupUtil.firstMarkupPoint;
                        h      = plot([origin(1);mouseMotionPoint(1)],[origin(2);mouseMotionPoint(2)],...
                            '-','Color',plotMarkupUtil.markupColour,'LineWidth',plotMarkupUtil.markupLineWidth);
                        
                        plotMarkupUtil.markupHandles    = [plotMarkupUtil.markupHandles; h];
                        plotMarkupUtil.firstMarkupPoint = [];
                    end
                end  
            end
        end

        % Callback for the mouse button up event.
        function mouseUpCallback(plotMarkupUtil,~,~)

            if plotMarkupUtil.isButtonDown
                
                plotMarkupUtil.isButtonDown = false;
                
                % Broadcast notice of event
                %notify(plotMarkupUtil,'PAUMouseButtonUp',IRUMouseButtonUpEvent(plotMarkupUtil.imgRoi)); 
                
            end  
        end        
    end
    
    methods (Access=private)
        
        % Displays the image in a figure.
        %
        % @param   plotMarkupUtil   The plot markup utility to show the
        %                           figure for.
        %
        function [plotMarkupUtil] = show(plotMarkupUtil)
            
            % create the figure
            plotMarkupUtil.plotFigHandle  = figure;
            plotMarkupUtil.plotAxisHandle = axes;

            % set the callback functions for the figure
            set(plotMarkupUtil.plotFigHandle, ...
                'WindowButtonDownFcn',   @plotMarkupUtil.mouseDownCallback, ...
                'WindowButtonUpFcn',     @plotMarkupUtil.mouseUpCallback,   ...
                'WindowButtonMotionFcn', @plotMarkupUtil.mouseMotionCallback);    
        end
        
        % Initializes the data members of the plot markup utility.
        %
        % @param   plotMarkupUtil   The plot markup instance to initialize.
        %
        % @return  The plot markup utility instance.
        %       
        function [plotMarkupUtil] = init(plotMarkupUtil)
            plotMarkupUtil.isButtonDown     = false;
            plotMarkupUtil.plotFigHandle    = [];
            plotMarkupUtil.plotAxisHandle   = [];
            plotMarkupUtil.markupHandles    = [];
            plotMarkupUtil.firstMarkupPoint = [];
            plotMarkupUtil.markupColour     = [1 0 0];
            plotMarkupUtil.markupLineWidth      = 0.5;
        end
        
        % Closes the currently open figure.
        %
        % @param   plotMarkupUtil   The plot markup instance to initialize.
        %
        % @return  The plot markup utility instance.
        %       
        function [plotMarkupUtil] = close(plotMarkupUtil)
            if ~isempty(plotMarkupUtil.plotFigHandle)                
                
                for i=1:length(plotMarkupUtil.markupHandles)
                    try
                        delete(plotMarkupUtil.markupHandles(i));
                    catch
                        % intentionally left blank
                    end
                end
                
                plotMarkupUtil.markupHandles = [];
                
                try
                    close(plotMarkupUtil.plotFigHandle);
                    plotMarkupUtil.plotFigHandle  = [];
                    plotMarkupUtil.plotAxisHandle = [];
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
        function [plotMarkupUtil] = validateColour(plotMarkupUtil,colour)
            
            if ischar(colour)
                colour = lower(colour);
                switch colour
                    case {'y','yellow'}
                        plotMarkupUtil.markupColour = [1 1 0];
                    case {'m','magenta'}
                        plotMarkupUtil.markupColour = [1 0 1];
                    case {'c','cyan'}
                        plotMarkupUtil.markupColour = [0 1 1];
                    case {'r','red'}
                        plotMarkupUtil.markupColour = [1 0 0];
                    case {'g','green'}
                        plotMarkupUtil.markupColour = [0 1 0];
                    case {'b','blue'}
                        plotMarkupUtil.markupColour = [0 0 1];
                    case {'w','white'}
                        plotMarkupUtil.markupColour = [1 1 1];
                    case {'k','black'}
                        plotMarkupUtil.markupColour = [0 0 0];
                    otherwise
                        printf('PlotMarkupUtility - colour string is invalid');
                end
            else
                if length(colour)~=3
                    fprintf('PlotMarkupUtility - colour vector must have 3 elements');
                elseif ~all(colour<=1 & colour>=0.0)
                    fprintf('PlotMarkupUtility - colour vector elements must be between 0 and 1');
                else 
                    plotMarkupUtil.markupColour = double(colour);
                end
            end
        end
    end
end

