%**************************************************************************
% Executes the Shift Cipher, encrypting the plaintext using the selected
% key OR decrypting the cipher text using the same key:
%
%                       e(x) = x + b mod m
%                       d(y) = y - b mod m
%
% @param   x   The plaintext message to encrypt (if non-empty).
% @param   y   The ciphertext message to decrypt (if non-empty).
% @param   b   The key to encrypt or decrypt.
% @param   m   The size of the symbol space.
%
% @return  The encrypted or decrypted message.
%
% @throw   Error if any message symbol element does not belong to the group 
%          Zm.
% @throw   Error if the key does not belong to the group Zm.
% @throw   Error if the symbol space size m is undefined (necessary for
%          non-string plain or encrypted text messages.
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
function [message] = cphr_shift(x,y,b,m)

    encrypt = 1;
    
    if exist('m','var')
        m = abs(ceil(m(1))); 
    end

    if ~isempty(x)
        ctype = class(x);
        if ischar(x)
            x = uint16(x);
            m = 256;
        end
    elseif ~isempty(y)
        ctype = class(y);
        if ischar(y)
            y = uint16(y);
            m = 256;
        end
        encrypt = 0;      
    else
        error('cphr_shift: symbol space size m is undefined');
    end

    b = abs(ceil(b(1)));

    if encrypt

        if b>=m
            error(['cphr_shift: key (%d) is greater than the symbol ' ...
                   'space size (%d)'],b,m);
        end            

        if isempty(find(x<0,1)) && isempty(find(x>=m,1))
            message = mod(x+b,m);
        else
            error('cphr_shift: plaintext message includes invalid values');
        end
        
    else

        if b>=m
            error(['cphr_shift: key (%d) is greater than the symbol ' ...
                   'space size (%d)'],b,m);
        end         
        
        if isempty(find(y<0,1)) && isempty(find(y>=m,1))
            message = mod(y-b,m);
        else
            error('cphr_shift: ciphertext message includes invalid values');
        end        
        
    end
    
    message = cast(message,ctype);
    
end

% usage:

% ptxt = 'A long time ago, in a galaxy far, far away...';
% key  = 11;
% ctxt=cphr_shift(ptxt,[],key)
% txt=cphr_shift([],ctxt,key)

