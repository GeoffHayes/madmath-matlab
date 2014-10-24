%**************************************************************************
% Renames all files in a given directory.
%
% @param   directory   The directory to rename all files to.
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
function rename_files(directory)

    % is the directory valid?
    if isdir(directory)
        
        % grab all files from this directory        
        fileList = dir(directory);
        
        % skip the first two "files" which just correpsond to '.' and '..'
        for i=4:size(fileList)
            
            fileSrc    = [directory '/' fileList(i).name];
            fileNoExt  = fileList(i).name(1:end-4);
            fileRepl   = ['ca' num2str(i-3)];
            fileDest   = strrep(fileSrc, fileNoExt, fileRepl);
            movefile(fileSrc, fileDest);
            
        end
    end
end

