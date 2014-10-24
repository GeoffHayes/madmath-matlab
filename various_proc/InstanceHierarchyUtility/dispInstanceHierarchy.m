function dispInstanceHierarchy(instance,level,name,prfx)

    if ~exist('prfx','var')
        prfx = '';
    end
    
    if ~exist('name','var')
        name = inputname(1);
    end

    if level>0
       
        fprintf([prfx '%s (%s)\n'],name,class(instance));
        
        if isstruct(instance) || isobject(instance)
            % get all fields from the struct
            fieldNames = fieldnames(instance);

            for i=1:length(fieldNames)
                dispInstanceHierarchy(instance.(fieldNames{i}),level-1,fieldNames{i},[prfx '\t']);
            end
        end 
    end
end