%**************************************************************************
% Creates the DES key schedule (for each iteration/round) given the 64 bit key.
%
% @param   key64     The 64-bit key.
% @param   encrypt   Flag indicating whether DES is encrypting or
%                    decrypting a message, as this affects the key
%                    schedule.
%
% @return  The key schedule.
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
function [keysched48] = des_createkeysched(key64,encrypt)

    % load the DES constants
    des_constants;
    
    global PC1;
    global PC2;
    global NUM_ROUNDS;
    
    % apply the permutation choice 1 to the 64-bit key to get a 56-bit key
    key56 = permutebits_u64(key64,PC1);
    
    % break the 56-bit key into halves, c and d of 28 bits each
    mask28 = cast(hex2dec('FFFFFFF'),'uint64');
    maskc  = cast(hex2dec('FFFFFFF0000000'),'uint64');
    maskd  = cast(hex2dec('0000000FFFFFFF'),'uint64');
    ci28   = bitshift(bitand(key56,maskc),-28);
    di28   = bitand(key56,maskd);
    
    keysched48 = cast(zeros(NUM_ROUNDS,1),'uint64');
    
    for i=1:NUM_ROUNDS
        switch i
            
            case {1,2,9,16}
                
                % shift bits of ci28 and di28 one position to left
                ci28 = bitand(bitor(bitshift(ci28,1),bitshift(bitand(ci28,hex2dec('8000000')),-27)),mask28);
                di28 = bitand(bitor(bitshift(di28,1),bitshift(bitand(di28,hex2dec('8000000')),-27)),mask28);
                
            otherwise
                
                % shift bits of ci28 and di28 two positions to left
                ci28 = bitand(bitor(bitshift(ci28,2),bitshift(bitand(ci28,hex2dec('C000000')),-26)),mask28);
                di28 = bitand(bitor(bitshift(di28,2),bitshift(bitand(di28,hex2dec('C000000')),-26)),mask28);
                
        end
        
        % combine ci28 and di28
        cidi56 = bitor(bitshift(ci28,28),di28);
        
        % apply the  permutation choice 2 to the 56-bit to get the 48-bit
        % key for round i (reverse the key schedule for decryption)
        if encrypt
            keysched48(i) = permutebits_u64(cidi56,PC2,56);
        else
            keysched48(NUM_ROUNDS-i+1) = permutebits_u64(cidi56,PC2,56);
        end
    end

end


