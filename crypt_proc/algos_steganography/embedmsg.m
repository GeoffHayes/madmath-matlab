%**************************************************************************
% Embeds a message within the least significant bits of a source.
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
function [target] = embedmsg(source, msg, numbits2emb)

    if ~strcmpi(class(msg),'uint8')
        error('embedmsg: message is not of type 8-bit unsigned integer');
    end
    
    if ~strcmpi(class(source),'uint8')
        error('embedmsg: source is not of type 8-bit unsigned integer');
    end
        
    if ndims(msg)>3
        error('embedmsg: message can only be at most three dimensional');
    end
    
    if ndims(msg)==3 && (size(msg,3)~=3 && size(msg,3)~=1)
        error('embedmsg: third dimension of message must be 1 or 3');
    end    
    
    if ~exist('numbits2emb','var')
        numbits2emb = 1;
    else
        numbits2emb = max(numbits2emb,4);
    end
    
    target = source;
    
    % determine the number of elements in the source (where each element
    % will have its LSB modified with a single bit from the message to
    % embed)
    numElemsSrc = numel(source);
    
    % note that numElemsSrc/8 is the number of bytes from the message that
    % we can embed 
    maxMsgBytesToEmbed = floor(numElemsSrc*numbits2emb/8);
    
    numBytesMsg = numel(msg);
    
    % determine the dimensions of the message
    [msgDims] = uint16(size(msg));
    
    while length(msgDims)<3
        msgDims = [msgDims 1];
    end
    
    % embed the message dimensions into the target
    atTargByte = 1;
    for u=1:3
        mask16 = uint16(1);
        for v=1:16
            if bitand(msgDims(u),mask16)
                target(atTargByte) = bitor(target(atTargByte),1);
            else
                target(atTargByte) = bitand(target(atTargByte),0);
            end
            
            mask16 = bitshift(mask16,1);
            
            atTargByte = atTargByte + 1;
        end
    end
    
    % reshape the message
    if ~isvector(msg)
       if ismatrix(msg)
           msg = reshape(msg',numBytesMsg,1);
       else
           % three dimensional
           d1 = msg(:,:,1)';
           d2 = msg(:,:,2)';
           d3 = msg(:,:,3)';

           msg = [d1(:) d2(:) d3(:)];
           
           msg = reshape(msg',numBytesMsg,1);
       end
    end

    % now embed as much of the message as we can into the target
    numBytesMsg = min(numBytesMsg,maxMsgBytesToEmbed-atTargByte+1);
    
    switch numbits2emb
        case 1
            mask = uint8(hex2dec('FE'));
        case 2
            mask = uint8(hex2dec('FC'));
        otherwise
            mask = uint8(hex2dec('F0'));
    end
    
    maskcmp = bitcmp(mask);

    for u=1:numBytesMsg
        
        msgByte = msg(u);

        % embed each bit of the message byte into a target byte
        for v=1:numbits2emb:8

            target(atTargByte) = bitor(bitand(target(atTargByte),mask),...
                                       bitand(msgByte,maskcmp));
                                   
            msgByte = bitshift(msgByte,-numbits2emb);                       
 
            atTargByte = atTargByte + 1;
        end
    end
end