%**************************************************************************
% Reads blocks of data of a known fixed size (nxn) with each having a text
% header.  Text file must be of format TEXT STRING - nxn block of data -
% TEXT STRING - nxn block of data - etc.
%
% @param   filename   The text file to read the data from.
% @param   n          The square block size.
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
function readblocks(filename,n)

    fid = fopen(filename, 'rt');

    % if file descriptor valid
    if fid>0
        
        % continue until end-of-file reached
        while ~feof(fid)
            % read/get the line of text (and ignore it)
            fgetl(fid);
            
            % read the nxn matrix of floating point data, transposing the 
            % result of the fscanf (since it populates A in column order)
            [A]    = fscanf(fid,'%f\n',[n,n])';
            
            % do something with the A
        end
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

