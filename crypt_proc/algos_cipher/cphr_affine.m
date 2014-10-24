%**************************************************************************
% Executes the Affine Cipher, encrypting the plaintext using the selected
% key OR decrypting the cipher text using the same key.
%
%                       e(x) = a*x + b mod m
%                       d(y) = a^-1*(y - b) mod m
%
% @param   x   The plaintext message to encrypt (if non-empty).
% @param   y   The ciphertext message to decrypt (if non-empty).
% @param   k   The key with components a and b to encrypt or decrypt.
% @param   m   The size of the symbol space.
%
% @return  The encrypted or decrypted message.
%
% @throw   Error if any message symbol element does not belong to the group 
%          Zm.
% @throw   Error if either key element (a,b) does not belong to the group Zm.
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
function [message] = cphr_affine(x,y,k,m)

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
        error('cphr_affine: symbol space size m is undefined');
    end
    
    a = abs(round(k.a(1)));
    b = abs(round(k.b(1)));

    if encrypt

        if a>=m || b>=m
            error(['cphr_affine: either a (%d) or b (%d) is greater than ' ...
                   'the symbol space size (%d)'],a,b,m);
        end            

        if isempty(find(x<0,1)) && isempty(find(x>=m,1))
            message = mod(a*x+b,m);
        else
            error('cphr_affine: plaintext message includes invalid values');
        end
        
    else

        if a>=m || b>=m
            error(['cphr_affine: either a (%d) or b (%d) is greater than ' ...
                   'the symbol space size (%d)'],a,b,m);
        end          
        
        if isempty(find(y<0,1)) && isempty(find(y>=m,1))
            ia      = ext_euclidean(a,m);
            message = mod(ia*(y-b),m);
        else
            error('cphr_affine: ciphertext message includes invalid values');
        end        
        
    end
    
    message = cast(message,ctype);
    
end

% usage:

% ptxt = 'A long time ago, in a galaxy far, far away...';
% key.a = 7;key.b = 3;
% ctxt=cphr_affine(ptxt,[],key)
% txt=cphr_affine([],ctxt,key)

