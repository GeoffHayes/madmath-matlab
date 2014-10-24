%**************************************************************************
% The sieve factoring algorithm returns all factors of the prime (if they
% exist) along with their exponents.
%
% @param   n    The integer to find a factors of.
%
% @return  A mx2 matrix of m factors (in column 1) along with their
%          corresponding exponents (in column 2).
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
function [F] = fctr_sieve(n)

    n  = uint64(round(abs(n(1))));
    
    F   = [];
    atf = 0;

    % use bit manipulation for even numbers
    mask      = uint64(1);
    num2fctrs = 0;
    while bitand(n,mask)==0
        num2fctrs = num2fctrs + 1;
        n         = bitshift(n,-1);
    end
    
    if num2fctrs>0
        atf      = atf+1;
        F(1,atf) = 2;
        F(2,atf) = num2fctrs;
    end
    
    % use bit manipulation for multiples of five
    mask      = uint64(5);
    num5fctrs = 0;
    while bitand(n,mask)==mask
        num5fctrs = num5fctrs + 1;
        n         = n/mask;
    end
    
    if num5fctrs>0
        atf      = atf+1;
        F(1,atf) = 5;
        F(2,atf) = num5fctrs;
    end 
    
    % use the sieve to get all primes from 2 to the (maybe) reduced n
    prms = sieve_erat(ceil(sqrt(double(n))));
    
    % remove the 2 and 5 primes since already accounted for
    prms(prms==2 | prms==5) = [];
    
    for j=length(prms):-1:1
        p         = prms(j);
        numpfctrs = 0;
        while ~mod(n,p)
            numpfctrs = numpfctrs+1;
            n         = n/p;
        end
        
        if numpfctrs>0
            atf      = atf+1;
            F(1,atf) = p;
            F(2,atf) = numpfctrs;
        end
    end
    
    if n~=1
        atf      = atf+1;
        F(1,atf) = n;
        F(2,atf) = 1;
    end
    
    F = sortrows(F',1);
   
end

