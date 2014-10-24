function [levelsDeep] = drawInstanceHierarchy(instance,level,name,h,x,y,xp,yp)

    levelsDeep = 1;

    if ~exist('h','var')
        h = figure;
        axis([0 100 0 100]);
        hold on;
        x = 10;
        y = 90;
    end
    
    if ~exist('xp','var')
        xp = [];
        yp = [];
    end
    
    if ~exist('name','var')
        name = inputname(1);
    end

    if level>0
       
        h = text(x,y,[name ' (' class(instance) ')'],'EdgeColor','b');
        ext = get(h,'Extent');
        
        % draw a line?
        if ~isempty(xp)
            line([xp xp],[yp ext(2)+ext(4)/2]);
            line([xp ext(1)],[ext(2)+ext(4)/2 ext(2)+ext(4)/2]);
        end
        
        x = ext(1)+ext(3)+2;
        y = ext(2)-5;
        levelsDeep = y;

        if isstruct(instance) || isobject(instance)
            % get all fields from the struct
            fieldNames = fields(instance);

            for i=1:length(fieldNames)

                levelsDeep = ...
                    drawInstanceHierarchy(instance.(fieldNames{i}),level-1,fieldNames{i},h,x,levelsDeep,ext(1)+ext(3)/2,ext(2));
            end
        end 
    end
end