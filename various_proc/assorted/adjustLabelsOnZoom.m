function adjustLabelsOnZoom
    close all;
    
    % pre-zoom x-axis limits
    PREZOOM_AXIS_X_LIMITS = [];
    % pre-zoom y-axis limits    
    PREZOOM_AXIS_Y_LIMITS = [];  

    % x and y coordinates of the labels
    X_COORDS=[10000;11000;12000];  
    Y_COORDS=[20000;21000;22000];  
    
    % labels
    points_id={'p100';'p200';'p300'};
    
    % scatter the points on the plot
    figure(1);
    scatter(X_COORDS, Y_COORDS, 'b^');
    grid on;
    
    % save handles to text/ids for post-zoom callback
    GRID_TEXT_HANDLES = text(X_COORDS, Y_COORDS+10, points_id); 
    
    % get the figure's zoom mode object
    h = zoom;
    
    % set the zoom callbacks
    set(h,'ActionPreCallback',@prezoom,'ActionPostCallback',@postzoom);
    
    function prezoom(~,event_obj)
        PREZOOM_AXIS_X_LIMITS = get(event_obj.Axes,'XLim');
        PREZOOM_AXIS_Y_LIMITS = get(event_obj.Axes,'YLim');
    end

    function postzoom(~,event_obj)

         % get the post-zoom limits for each axis
         PSTZOOM_AXIS_X_LIMITS = get(event_obj.Axes,'XLim');
         PSTZOOM_AXIS_Y_LIMITS = get(event_obj.Axes,'YLim');
         
         % determine the difference in the limits
         prezoomYDiff = abs(PREZOOM_AXIS_Y_LIMITS(2) - PREZOOM_AXIS_Y_LIMITS(1));
         pstzoomYDiff = abs(PSTZOOM_AXIS_Y_LIMITS(2) - PSTZOOM_AXIS_Y_LIMITS(1));
         
         prezoomXDiff = abs(PREZOOM_AXIS_X_LIMITS(2) - PREZOOM_AXIS_X_LIMITS(1));
         pstzoomXDiff = abs(PSTZOOM_AXIS_X_LIMITS(2) - PSTZOOM_AXIS_X_LIMITS(1));
         
         % adjust each label
         for i=1:length(Y_COORDS)
            % get the coordinate only for the marker
            y = Y_COORDS(i);
            x = X_COORDS(i);
            
            % get the text handle for the label
            txtHandle = GRID_TEXT_HANDLES(i);
            
            % get the 3D position of the label
            txtPos    = get(txtHandle,'Position');
            
            % calculate the difference (between the marker and label) in x
            % and y
            yDiff     = abs(txtPos(2)-y);
            xDiff     = abs(txtPos(1)-x);
            
            % set the new difference/range (from marker to label) to be at least 1
            % and at most 10 (chosen as per original placement) using the ratio of
            % post-zoom y-axis limit wrt pre-zoom y-axis limit
            newYDiff  = min(10,max(1,pstzoomYDiff/prezoomYDiff*yDiff));
            newXDiff  = min(10,max(0,pstzoomXDiff/prezoomXDiff*xDiff));
            
            
            % set the new position for the label relative to the marker
            set(txtHandle,'Position', [x+newXDiff y+newYDiff txtPos(3)]);
         end
    end
end