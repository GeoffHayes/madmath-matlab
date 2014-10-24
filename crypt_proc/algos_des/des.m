%**************************************************************************
% Executes the 64-bit Data Encryption Standard (DES) algorithm.
%
% @param   xstr   The plaintext message (stream) to encrypt (if non-empty).
% @param   ystr   The ciphertext message (stream) to decrypt (if non-empty).
% @param   k      The 64-bit key used for encryption and decryption.
%
% @return  The encrypted or decrypted message.
%
% @throw   Error if the key is not a hexadecimal string.
% @throw   Error if the data type is unsupported.
%
% @warning The key is padded with zeros if the input key does not contain
%          16 characters.  The key is truncated if the input key contains
%          more than 16 characters.
%
% [1] Stinson, Douglas R., Cryptography Theory and Practice, CRC Press
%     1995.
% [2] http://en.wikipedia.org/wiki/DES_supplementary_material
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
function [message] = des(xstr,ystr,k)

    % load the DES constants
    des_constants;
    
    global IP; 
    global IPi; 
    global NUM_BITS_PER_BYTE;
    global NUM_BYTES_PER_BLK;
    global MASK_LEFT;
    global MASK_RIGHT;
    global HALF_BLK_SIZE;
    global NUM_ROUNDS;

    % define the key up to 64 bits
    key64 = uint64(0);

    if ischar(k)
        try
            % key string is marked as hexadecimal
            strlen = min(16,length(k));
            for i=1:strlen
                key64 = bitshift(key64,4);
                key64 = bitor(key64,hex2dec(k(i)));
            end
        catch
            error('des: key contains character that is invalid (%s)',k(i));
        end
    else
        error('des: key is not a hexidecimal string');
    end

    % read all data from the stream and encrypt or decrypt as appropriate
    if ~isempty(xstr)
       encrypt = 1;
       istr    = xstr;
    else
       encrypt = 0;
       istr    = ystr;
    end
    
    % create the key schedule
    keysched = des_createkeysched(key64,encrypt);
    
    ctype = class(istr);

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
    elseif ischar(istr)
        numbits = 8;
    else
        error('des: unsupported class type %s',ctype);
    end   
    
    numbytes      = numbits/NUM_BITS_PER_BYTE;
    numelems      = numel(istr);
    numinblks2add = NUM_BYTES_PER_BLK/numbytes;

    % may need to pad the output message array as DES encrypts every 64
    % bits worth of data
    xtra = rem(numelems,numinblks2add);
    if xtra~=0
        data2pad = numinblks2add-xtra;
    else
        data2pad = 0;
    end
    
    message       = cast(zeros(1,numelems+data2pad),ctype);
    
    OUTBLK_MASK   = cast(hex2dec([repmat('F',1,numbits/4) repmat('0',1,16-numbits/4)]),'uint64');
    
    % convert the input stream data to unsigned 64 bits to avoid future
    % casts
    istr = uint64(istr);
    
    i=1;

    while i<=numelems
        
        % initialize a block of data
        block = uint64(0);
        
        bytesadded = 0;
        
        % read enough in-stream data to fill the block
        strtbyte = i;
        stopbyte = min(numelems,numinblks2add+i-1);
        for i=strtbyte:stopbyte
            block      = bitor(bitshift(block,numbits),istr(i));
            bytesadded = bytesadded + 1;
        end
        i=i+1;

        % pad with zeros if not enough bytes have been added
        block = bitshift(block,(numinblks2add-bytesadded)*NUM_BITS_PER_BYTE);
        
        % permute the block according to the initial permutation matrix IP
        block = permutebits_u64(block,IP);
        
        % divide the block into its left and right halves ('p'=='previous')
        ljp32 = bitshift(bitand(block,MASK_LEFT),-HALF_BLK_SIZE);
        rjp32 = bitand(block,MASK_RIGHT);
        
        % do the iterations
        for j=1:NUM_ROUNDS
            lj32  = rjp32;
            rj32  = bitxor(ljp32,applyf(rjp32,keysched(j)));
            rjp32 = rj32;
            ljp32 = lj32;
        end
    
        % combine the left and right halves
        block = bitor(bitshift(rjp32,HALF_BLK_SIZE),ljp32);
        
        % permute the block according to the inverse permutation matrix IPi
        block = permutebits_u64(block,IPi);
        
        % write out the block
        numbytes2ext = numinblks2add;
        mask         = OUTBLK_MASK;
        
        for k=strtbyte:numinblks2add+strtbyte-1
            message(k)   = bitshift(bitand(block,mask),-(numbytes2ext-1)*numbits);
            mask         = bitshift(mask,-numbits);
            numbytes2ext = numbytes2ext -1;
        end
        
    end

end

%**************************************************************************
% Executes the DES f-function against the 32-bit block with a 48-bit key.
%
% @param   iblk32   Theinput 32-bit integer.
% @param   key48    The 48-bit key.
%
% @return  The output 32-bit block.
%
% [1] Stinson, Douglas R., Cryptography Theory and Practice, CRC Press
%     1995.
% [2] http://en.wikipedia.org/wiki/DES_supplementary_material
%**************************************************************************
function [oblk32] = applyf(iblk32,key48)

    global E;  
    global MASK_B1;
    global MASK_B2345;
    global MASK_B6;
    global SBoxes;
    global P;
    
    % expand the input block to 48-bits
    iblk48 = permutebits_u64(iblk32,E,32);
    
    % exclusive or with the key
    iblk48 = bitxor(iblk48,key48);

    % for each 8-bit block extract the six bits and replace with that from
    % the ith S-box
    maskb1    = MASK_B1;
    maskb2345 = MASK_B2345;
    maskb6    = MASK_B6;
    
    oblk32 = cast(0,'uint64');
    
    % define the shift (to the left)
    shft = -6;
    
    for i=1:8
        
        % extract the bits for the row
        row = bitshift(bitand(iblk48,maskb1),shft*(8-i)-4);
        row = bitor(row,bitshift(bitand(iblk48,maskb6),shft*(8-i))) + 1;
        
        % extract the bits for the column
        col = bitshift(bitand(iblk48,maskb2345),shft*(8-i)-1) + 1;
        
        % get the entry from the ith S-box, shifting the output block
        % appropriate number of bits to the right (a function of 4-bits
        % since the output from the SBox is a 4-bit integer)
        oblk32  = bitor(bitshift(oblk32,4),SBoxes(i,row,col));
        
        % shift the masks so that we are ready to extract the next set of
        % 6-bits
        maskb1    = bitshift(maskb1,shft);
        maskb2345 = bitshift(maskb2345,shft);
        maskb6    = bitshift(maskb6,shft);
        
    end
    
    % apply the fixed permutation P
    oblk32 = permutebits_u64(oblk32,P,32);

end

% usage:

% ptxt='A long time ago in a galaxy far, far away......';
% ctxt=des(ptxt,[],'133457799bbcdff1');
% ictxt=des([],ctxt,'133457799bbcdff1');
