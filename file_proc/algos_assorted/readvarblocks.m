%**************************************************************************
% Reads blocks of data of a variable size with each having a text header.  
% Text file must be of format TEXT STRING - mxn block of data -
% TEXT STRING - uxv block of data - etc.
%
% @param   filename   The text file to read the data from.
% @param   tag        The tag string that is common to all text headers.
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
function readvarblocks(filename,tag)

    fid = fopen(filename, 'rt');

    % if file descriptor valid
    if fid>0
        
        data = [];

        % continue until end-of-file reached
        while ~feof(fid)

            % read/get the current line of text
            currentline = fgetl(fid);

            if ~isempty(strfind(currentline,tag))
                
                if ~isempty(data)
                   % do something with the mxn block of data 
                end
                
                % since current line has the tag, then move to the next line
                currentline = fgetl(fid);
                
                data = [];
                at   = 1;
            end
            % since pattern is line of text then numbers, then can assume
            % that the current line is a set of numbers separated by spaces
            % so convert to a numeric array
            data(at,:) = str2num(currentline);
            at         = at + 1;
            
        end
        
        if ~isempty(data)
        	% do something with the mxn block of data 
        end
                

        % close the file
        fclose(fid);
    end
end

% sample data
% # TT50/Data/BaseballPitch/v_BaseballPitch_g01_c01
%  5 0.25 0.228125 0.0654206
%  5 0.133333 0.0375 0.0747664
%  5 0.208333 0.55625 0.0747664
%  5 0.495833 0.221875 0.0747664
%  #TT50/Data/BaseballPitch/v_BaseballPitch_g01_c02
%  5 0.591667 0.134375 0.0860215
%  5 0.320833 0.125 0.0967742
%  5 0.458333 0.24375 0.0967742
%  5 0.520833 0.140625 0.0967742
%  # TT50/Data/BaseballPitch/v_BaseballPitch_g01_c03
%  5 0.625 0.821875 0.0873786
%  5 0.6125 0.765625 0.203883
%  5 0.575 0.78125 0.262136
%  5 0.6 0.778125 0.271845

