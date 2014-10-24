%**************************************************************************
% Executes the RSA decryption algorithm against a stream of data for a
% private key.
%
% @param   ystr    The ciphertext message (stream) to decrypt.
% @param   prikey  The private RSA key.
%
% @return  The descrypted message.
%
% @throw   Error if the data type is unsupported.
%
% [1] Stinson, Douglas R., Cryptography Theory and Practice, CRC Press
%     1995.
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
function [xstr] = rsa_dec(ystr,prikey)

    error('rsa_dec: incomplete, need to fix the reading of 16bit blocks and writing of 8bit blocks');

    % load the RSA constants
    rsa_constants;
    
    global NUM_BITS_PER_BYTE;
    global NUM_BYTES_PER_EBLK;
    
    ctype = class(ystr);

    if strcmpi(ctype,'uint64')     || strcmpi(ctype,'int64') 
        numbits     = 64;
    elseif strcmpi(ctype,'uint32') || strcmpi(ctype,'int32')
        numbits = 32;
    elseif strcmpi(ctype,'uint16') || strcmpi(ctype,'int16')
        numbits     = 16;
    elseif strcmpi(ctype,'uint8')  || strcmpi(ctype,'int8')
        numbits     = 8;
    elseif strcmpi(ctype,'double')
        numbits     = 64;
    elseif strcmpi(ctype,'single')
        numbits = 32;
    elseif ischar(ystr)
        numbits = 8;
    else
        error('rsa_dec: unsupported class type %s',ctype);
    end   
    
    numbytes      = numbits/NUM_BITS_PER_BYTE;
    numelems      = numel(ystr);
    numinblks2add = NUM_BYTES_PER_EBLK/numbytes;

    % may need to pad the output message array
    xtra = rem(numelems,numinblks2add);
    if xtra~=0
        data2pad = numinblks2add-xtra;
    else
        data2pad = 0;
    end
    
    xstr = cast(zeros(1,(numelems+data2pad)/2),ctype);
    
    OUTBLK_MASK = cast(hex2dec([repmat('F',1,numbits/4) ...
        repmat('0',1,16-numbits/4)]),'uint64');
    
    % convert the input stream data to unsigned 64 bits to avoid future
    % casts
    ystr = uint64(ystr);
    
    i=1;
    
    while i<=numelems
        
        % reset the block
        block64 = uint64(0);

        bytesadded = 0;
        
        % read enough in-stream data to fill the block
        strtbyte = i;
        stopbyte = min(numelems,numinblks2add+i-1);
        for i=strtbyte:stopbyte
            block64 = bitor(bitshift(block64,numbits),...
                ystr(i));
            bytesadded = bytesadded + 1;
        end
        i=i+1;

        % pad with zeros if not enough bytes have been added
        block64 = bitshift(block64,...
            (numinblks2add-bytesadded)*NUM_BITS_PER_BYTE);
        
        % encrypt the block
        block64 = safe_moduloexp(block64,prikey.a,prikey.n);

        % write out the block
        numbytes2ext = numinblks2add;
        % since we only write out the last 32 bits, shift the mask
        mask         = bitshift(OUTBLK_MASK,-48);
        
        for k=strtbyte:numinblks2add+strtbyte-1
            xstr(k)      = bitshift(bitand(block64,mask),...
                -(numbytes2ext-1)*numbits);
            mask         = bitshift(mask,-numbits);
            numbytes2ext = numbytes2ext -1;
        end
    end
end



% usage:

% ptxt='A long time ago in a galaxy far, far away......';
% [pub,pri]=rsa_generatekeys(9,1003)
% ctxt=rsa_enc(ptxt,pub);
