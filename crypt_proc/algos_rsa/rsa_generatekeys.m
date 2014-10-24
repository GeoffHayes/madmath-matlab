%**************************************************************************
% Generates the private and public keys for the RSA cyrptosystem given an
% integer range to draw the (secret) prime numbers (p,q) that are used to
% determine n (=p*q).
%
% @param   lb   The lower bound on the allowed prime numbers.
% @param   ub   The upper bound on the allowed prime numbers.
%
% @return  The public and private keys.
%
% @throw   Error if no two prime numbers can be found given the lower and
%          upper bounds.
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
function [pubkey,prikey] = rsa_generatekeys(lb,ub)

    rsa_constants;
    
    global MAX_PQ;
    global MIN_PRIME_NUMBER;

    % load the prime numbers
    load prime_numbers.mat;

    MAX_ITERS  = 1000;
    
    lb = round(abs(lb(1)));
    ub = round(abs(ub(1)));
    
    if lb > ub
        tmp = ub;
        ub  = lb;
        lb  = tmp;
    end
    
    lb = max(lb,MIN_PRIME_NUMBER);
    ub = max(lb,ub);
 
    len   = length(prime_numbers);
    lbidx = [];
    ubidx = [];

    for i=1:len
        if prime_numbers(i)>=lb
            lbidx=i;
            break;
        end
    end
    
    if isempty(lbidx)
        error('rsa_generatekeys: no prime number at least %d found in db',lb);
    end
    
    for i=lbidx:len
        if prime_numbers(i)>=ub
            if prime_numbers(i)==ub
                ubidx=i;
            elseif i>1
                ubidx=i-1;
            end
            break;
        end
    end
    
    if isempty(ubidx)
        error('rsa_generatekeys: no prime number at most %d found in db',ub);
    end
    
    if ubidx<lbidx
        error('rsa_generatekeys: no prime numbers within the interval [%d,%d]',lb,ub);
    end
    
    % generate two prime numbers given the indices into the interval
    if (ubidx-lbidx)==0
        % only one possible prime number for both p and q
        p = prime_numbers(lbidx);
        q = p;
    else
        p = prime_numbers(lbidx+round((ubidx-lbidx)*rand));
        q = p;
        while q==p
            q = prime_numbers(lbidx+round((ubidx-lbidx)*rand));
        end
    end
    
    % ensure that p and q are such that their product does not exceed the
    % maximum allowed for a 64-bit unsigned integer
    iter = 0;
    while p>MAX_PQ/q
        iter = iter + 1;
        if (ubidx-lbidx)==0
            % only one possible prime number for both p and q
            error('rsa_generatekeys: could not generate valid p,q pair');;
        else
            p = prime_numbers(lbidx+round((ubidx-lbidx)*rand));
            q = p;
            while q==p
                q = prime_numbers(lbidx+round((ubidx-lbidx)*rand));
            end
        end
        if iter==MAX_ITERS
            error('rsa_generatekeys: could not generate valid p,q pair');
        end
    end
        
    pubkey.n = p*q;
    prikey.p = p;
    prikey.q = q;
    prikey.n = pubkey.n;
    
    phin = (p-1)*(q-1);
    
    % generate a random number b in the interval 1...phin-1
    pubkey.b = 1 + round((phin-1)*rand);
    while euclidean(pubkey.b,phin)~=1
        pubkey.b = 1 + round((phin-1)*rand);
    end
    
    % generate the inverse of b
    prikey.a = ext_euclidean(pubkey.b,phin);
    
    clear prime_numbers;

end

