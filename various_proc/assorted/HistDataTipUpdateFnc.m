function output_txt = HistDataTipUpdateFnc(~,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

 % get the position data
 pos = get(event_obj,'Position');
 
 % start the tip text
 output_txt = {['Bin Count: ',num2str(pos(2),4)],''};
 
 % get the xdata to determine the bar stats
 xdata = get(event_obj.Target,'XData');
 
 % get the child handles to see how many bars there are per bin (need
 % to handle multi-bar per bin differently the single bar per bin case)
 childhandles=get(get(event_obj.Target,'Parent'),'Children');
 
 if ~isempty(xdata) && ~isempty(childhandles)
    % get the bar width, min and max values
    minx  = min(xdata(:));
    maxx  = max(xdata(:));
    diff  = xdata(1,2)-xdata(1,1);
    
    numchildren = length(childhandles);
    
    lowbnd = 0;
    uppbnd = 0;
    binctr = 0;
    
    if numchildren==1
        % only one bar per bin, so easy to determine the lower and upper bounds for
        % the bar
        lowbnd = pos(1)-diff/2;
        uppbnd = pos(1)+diff/2;
        binctr = pos(1);
    else
        % else more than one bar per bin
        % note that left-most bar corresponds to the last child, and
        % right-most bar corresponds to the first child
        % get the far edges of the last bar in each set of numchildren bars
        faredges = get(childhandles(1),'XData');faredges=faredges(3,:);
        
        % get the near edges of the first bar in each set of numchildren
        % bars
        nearedges = get(childhandles(end),'XData');nearedges=nearedges(1,:);
        
        % determine the centres for each group of numchildren bars
        barcentres = (faredges+nearedges)/2;
        
        % determine the intervals for each bar centre
        if length(barcentres)>1
            barintervals = [-Inf (barcentres(2:end)+barcentres(1:end-1))/2 Inf];
        elseif length(barcentres)==1
            barintervals = [-Inf,Inf];
        end
        
        % find the lower and upper bounds for the data
        ubidx = find(barintervals>=pos(1),1);
        lbidx = ubidx-1;
        
        lowbnd = barintervals(lbidx);
        uppbnd = barintervals(ubidx);
        
        binctr = barcentres(lbidx);
    end
    
    if lowbnd<=minx
        lowbnd = -Inf;
    end
    if uppbnd>=maxx
        uppbnd = Inf;
    end
    % finish the tip text
    output_txt = [output_txt {['Bin Center: ', num2str(binctr,4)], ...
                              ['Bin Edges: [', num2str(lowbnd,4), ...
                               ', ', num2str(uppbnd,4),']']}];
 end

