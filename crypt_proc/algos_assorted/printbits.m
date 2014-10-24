%**************************************************************************
% Prints the bits in the passed number.
%
% @param   x      The number to print the bits for.
% @param   skip   Flag indicating whether every skipth bit should be
%                 ignored.
%
% @throw   Error if the data type is unsupported.
% @throw   Error if x is not a single value.
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
function printbits(x,skip)

    if length(x)~=1
        error('printbits: only single value inputs are supported');
    end

    ctype = class(x);

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
    elseif ischar(x) 
        numbits = 8;
    else
        error('printbits: unsupported class type %s',ctype);
    end
    
    x    = cast(x,'uint64');
    mask = cast(2^(numbits-1),'uint64');
    
    if exist('skip','var')
        for i=1:numbits
            if mod(i,skip)>0
                if bitand(x,mask)
                    fprintf('1');
                else
                    fprintf('0');
                end
            end
            mask = bitshift(mask,-1);
        end
    else
        for i=1:numbits
            if bitand(x,mask)
                fprintf('1');
            else
                fprintf('0');
            end
            mask = bitshift(mask,-1);
        end
    end
    fprintf('\n');
end