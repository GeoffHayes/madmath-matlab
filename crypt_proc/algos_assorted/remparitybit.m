%**************************************************************************
% Removes the parity bit at the end of each byte.
%
% @param   xstr     The "stream" of bytes to remove the parity bit from.
%
% @note    This function will reduce the number of set bits within an
%          integer.  For example, a 64-bit unsigned int will have 8 parity
%          bits removed (as there are 8 bytes) leaving 56 set bits.
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
function[xstr] = remparitybit(xstr)

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
    numvals  = numel(xstr);
    numbytes = numbits/NUM_BITS_PER_BYTE;
    
    for i=1:numvals
        
        x    = xstr(i);
        y    = cast(0,ctype);
        mask = cast(bitshift(hex2dec('FFFE'),(numbytes-1)*NUM_BITS_PER_BYTE),ctype);
        
        for j=1:numbytes

            % extract the seven bits (of byte j)
            y    = bitor(y,bitshift(bitand(x,mask),-(numbytes-j+1)));
            mask = bitshift(mask,-NUM_BITS_PER_BYTE);

        end
        
        xstr(i) = y;
        
    end
end