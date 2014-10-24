%**************************************************************************
% Executes the Substution Cipher, encrypting the plaintext using the 
% selected key OR decrypting the cipher text using the same key.
%
%                       e(x) = p(x)
%                       d(y) = ip(y)
%
% @param   x   The plaintext message to encrypt (if non-empty).
% @param   y   The ciphertext message to decrypt (if non-empty).
% @param   p   The key to encrypt or decrypt (permutation of all possible 
%              symbols).
% @param   m   The size of the symbol space.
%
% @return  The encrypted or decrypted message.
%
% @throw   Error if any message symbol element does not belong to the group 
%          Zm.
% @throw   Error if the key (permutation) does not include all symbols in
%          the symbol space.
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
function [message] = cphr_subs(x,y,p,m)

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
        error('cphr_subs: symbol space size m is undefined');
    end

    if encrypt

        if length(p)~=m || length(unique(p))~=m || min(p)~=0 | max(p)~=(m-1)
            error('cphr_subs: key is not a unique permutation of all symbols');
        end            

        if isempty(find(x<0,1)) && isempty(find(x>=m,1))
            
            message = zeros(size(x));
            for i=1:length(message)
               message(i) = p(x(i)+1); 
            end
            
        else
            error('cphr_subs: plaintext message includes invalid values');
        end
        
    else

        if length(p)~=m || length(unique(p))~=m || min(p)~=0 | max(p)~=(m-1)
            error('cphr_subs: key is not a unique permutation of all symbols');
        end            
        
        if isempty(find(y<0,1)) && isempty(find(y>=m,1))
            
            % invert the key
            ip = zeros(size(p));
            for i=1:length(ip)
                ip(p(i)+1) = i-1;
            end

            message = zeros(size(y));
            for i=1:length(message)
               message(i) = ip(y(i)+1); 
            end            
            
        else
            error('cphr_subs: ciphertext message includes invalid values');
        end        
        
    end
    
    message = cast(message,ctype);
    
end

% usage:
% ptxt = 'A long time ago, in a galaxy far, far away...';
% key=randperm(255)-1;
% ctxt = cphr_subs(ptxt,[],key)
% text = cphr_subs([],ctxt,key)
