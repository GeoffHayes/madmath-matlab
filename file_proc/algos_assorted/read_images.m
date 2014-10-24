%**************************************************************************
% Reads all images from a directory into a cell array.  The cell array is
% then saved back to the directory as a images.mat file (overwriting any
% pre-existing file of the same name).
%
% @param   directory   The directory to read all the images from.
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
function read_images(directory)

    % verify that the directory is valid
    if isdir(directory)
        
        % grab all files from this directory        
        files = dir(directory);
        
        % find the number of files in the directory
        numimgs = length(files);
        
        % intialize the cell array of images
        images = cell(numimgs,1);
        
        % read all files
        k = 0;
        for i=1:numimgs
            
            % append the file name to the directory
            filesrc = fullfile(directory,files(i).name);
            
            % only read into the workspace those images that have a certain
            % extension
            [~, ~, ext] = fileparts(filename);
            
            switch ext
                case {'jpg','gif','bmp','tif','png'}
                    k=k+1;
                    images{k}=imread(filesrc);
                otherwise
                    % do nothing
            end
        end
        
        % remove those empty elements in the cell array (if any exist)
        if k~=numimgs
            images = images(1:k);
        end
        
        % save the cell array to the same directory
        save(fullfile(directory,'images.mat'),'images'); 
        
    end
end
