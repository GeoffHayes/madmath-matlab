%**************************************************************************
% Adds a parity bit at the end of each byte indicating the eveness or
% oddness of the number of ones in the byte.
%
% @param   xstr     The "stream" of bytes to add a parity bit to.
% @param   doeven   Flag indicating whether the parity bit should be one if
%                   the number of ones in the byte is odd, else should be
%                   zero if the number of ones in the byte is even.
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
function[xstr] = addparitybit(xstr,doeven)

    if ~exist('doeven','var')
        doeven = 1;
    end

    NUM_BITS_PER_BYTE = 8;

    ctype = class(xstr);

    if strcmpi(ctype,'uint64')     || strcmpi(ctype,'int64') 
        numbits = 64;
    elseif strcmpi(ctype,'uint32') || strcmpi(ctype,'int32')
        numbits = 32;
    elseif strcmpi(ctype,'uint16') || strcmpi(ctype,'int16')
        numbits = 16;
    elseif strcmpi(ctype,'uint8')  || strcmpi(ctype,'int8')
        numbits = 8;
    elseif strcmpi(ctype,'double')
        numbits = 64;
    elseif strcmpi(ctype,'single')
        numbits = 32;
    elseif ischar(xstr)
        numbits = 8;
    end
    
    % loop through all elements and add the parity bit at the end of each
    % bite
    numvals = numel(xstr);
    
    numbytes = numbits/NUM_BITS_PER_BYTE;
    
    for i=1:numvals
        
        y = xstr(i);
        
        mask = cast(bitshift(1,numbits-1),ctype);
        for j=1:numbytes

            numones = 0;
            
            % do all bits less one
            for k=1:NUM_BITS_PER_BYTE-1
                if bitand(mask,y)
                    numones = numones + 1;
                end
                mask = bitshift(mask,-1);
            end
            
            % do the final bit separately to delay the bitshift
            if bitand(mask,y)
            	numones = numones + 1;
            end

            if (doeven && mod(numones,2)==0) || (~doeven && mod(numones,2)~=0)
                % set the parity bit to be zero if not already zero
                if bitand(y,mask)~=0
                    y = bitxor(y,mask);
                end
            else
                % set the parity bit to be one if not already one
                if bitand(y,mask)==0
                    y = bitor(y,mask);
                end
            end
            
            mask = bitshift(mask,-1);

        end
        
        xstr(i) = y;
        
    end
end