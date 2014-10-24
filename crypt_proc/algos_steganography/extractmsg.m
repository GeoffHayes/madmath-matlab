%**************************************************************************
% Extracts a message from the least significant bits of a source.
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
function [msg] = extractmsg(source, numbits2ext)
 
    if numel(source)<16
        error('extactmsg: source does not have enough data');
    end
    
    if ~strcmpi(class(source),'uint8')
        error('extactmsg: source is not of type 8-bit unsigned integer');
    end    
    
    if ~exist('numbits2ext','var')
        numbits2ext = 1;
    else
        numbits2ext = max(numbits2ext,4);
    end
    
    % determine the number of elements in the source (where each element
    % has had its LSB modified with a single bit from the message to
    % extract)
    numElemsSrc = numel(source);

    % determine the dimensions of the message
    msgDims = uint16(zeros(1,3));
    
    % extract the message dimensions from the source
    atSrcByte = 1;
    for u=1:3
        mask16 = uint16(1);
        for v=1:16
            if bitand(source(atSrcByte),1)
                msgDims(u) = bitor(msgDims(u),mask16);
            end
            
            mask16 = bitshift(mask16,1);
            
            atSrcByte = atSrcByte + 1;
        end
    end
    
    numElemsMsg = prod(msgDims);
    
    % allocate memory to the message
    msg = uint8(zeros(numElemsMsg,1));
    
    switch numbits2ext
        case 1
            mask = uint8(hex2dec('FE'));
        case 2
            mask = uint8(hex2dec('FC'));
        otherwise
            mask = uint8(hex2dec('F0'));
    end
    
    maskcmp = bitcmp(mask);    
    
    for u=1:length(msg)

        for v=1:numbits2ext:8
            
            msg(u) = bitor(msg(u),bitshift(bitand(source(atSrcByte),maskcmp),(v-1)));

            atSrcByte = atSrcByte + 1;
            
            if atSrcByte>numElemsSrc
                break;
            end
        end
        
        if atSrcByte>numElemsSrc
            break;
        end
    end 
    
    % reshape the message
    if msgDims(3)==3
        d1 = reshape(msg(1:3:end),msgDims(2),msgDims(1));
        d2 = reshape(msg(2:3:end),msgDims(2),msgDims(1));
        d3 = reshape(msg(3:3:end),msgDims(2),msgDims(1));
        
        msg = uint8(zeros(msgDims));
        msg(:,:,1) = d1';
        msg(:,:,2) = d2';
        msg(:,:,3) = d3';
    else
        msg = reshape(msg,msgDims(2),msgDims(1))';
    end
end