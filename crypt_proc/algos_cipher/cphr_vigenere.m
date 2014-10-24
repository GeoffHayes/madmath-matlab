%**************************************************************************
% Executes the Vigenere Cipher, encrypting the plaintext using the selected
% key OR decrypting the cipher text using the same key.
%
%               e(x1,x2,..,xn) = (x1+k1,x2+k2,...,xn+kn) mod m
%               d(x1,x2,..,xn) = (x1-k1,x2-k2,...,xn-kn) mod m
%
% @param   x   The plaintext message to encrypt (if non-empty).
% @param   y   The ciphertext message to decrypt (if non-empty).
% @param   k   The key of length n.
% @param   m   The size of the symbol space.
%
% @return  The encrypted or decrypted message.
%
% @throw   Error if any message symbol element does not belong to the group 
%          Zm.
% @throw   Error if any key element does not belong to the group Zm.
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
function [message] = cphr_vigenere(x,y,k,m)

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
        error('cphr_vigenere: symbol space size m is undefined');
    end
    
    % convert the string key to numbers (if necessary)
    if ischar(k)
        k = double(k);
    end

    if ~isempty(find(k<0,1)) || ~isempty(find(k>=m,1))
        error('cphr_vigenere: key includes invalid values');
    end
    
    blocksize = length(k);
    
    if encrypt
        
        if isempty(find(x<0,1)) && isempty(find(x>=m,1))
            
            message = zeros(size(x));
            msglen  = length(x);
            
            iters = floor(msglen/blocksize);
            at     = 1;
            
            for i=1:iters
                message(at:at+blocksize-1) = mod(x(at:at+blocksize-1)+k,m);
                at = at+blocksize;
            end
            
            if at<=msglen
                message(at:end) = mod(x(at:end)+k(1:msglen-at+1),m);
            end
            
        else
            error('cphr_vigenere: plaintext message includes invalid values');
        end
        
    else         
        
        if isempty(find(y<0,1)) && isempty(find(y>=m,1))

            message = zeros(size(y));
            msglen  = length(y);
            
            iters = floor(msglen/blocksize);
            at     = 1;
            
            for i=1:iters
                message(at:at+blocksize-1) = mod(y(at:at+blocksize-1)-k,m);
                at = at+blocksize;
            end
            
            if at<=msglen
                message(at:end) = mod(y(at:end)-k(1:msglen-at+1),m);
            end            
            
        else
            error('cphr_vigenere: ciphertext message includes invalid values');
        end        
        
    end
    
    message = cast(message,ctype);
    
end

% usage:

% ptxt = 'A long time ago, in a galaxy far, far away...';
% key  = 'trek';
% ctxt=cphr_vigenere(ptxt,[],key)
% txt=cphr_vigenere([],ctxt,key)

