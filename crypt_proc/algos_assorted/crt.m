%**************************************************************************
% Executes the Chinese Remainder Theorem which determines a unique solution
% for the system of r congruences x=mod(ai,mi) for i=1,2,..,r where
% m1,m2,...,mr are pairwise relatively prime integers, and all ai are
% integers too.  The unique solution modulo M=m1*m2*...*mr is given by:
% 
%                   x = SUM(i=1..r)mod(ai*Mi*yi,M)
%
% @param   m   A vector of r integers as defined in the above equation.
% @param   a   A vector of r integers as defined in the above equation.
%
% @return  The unique solution modulo M.
%
% @throw   Error if m is not a vector.
% @throw   Error if a is not a vector.
% @throw   Error if m and a are not the same size.
% @throw   Error if all elements in m are not pairwise relatively prime.
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
function [x] = crt(m,a)

    if ~isvector(m)
        error('crt: m is not a vector');
    elseif ~isvector(a)
        error('crt: a is not a vector');
    elseif length(m)~=length(a)
        error('crt: m and a are not of the same length');
    end
    
    m = floor(m);
    a = floor(a);
    
    if ~isrelprime(m)
        error('crt: m contains at least two elements that are not relatively prime');
    end
    
    len  = length(m);
    M    = zeros(size(m));
    y    = M;
    Mprd = 1;

    for i=1:len
        Mprd = Mprd*m(i);
    end
    
    for i=1:len
        M(i) = Mprd/m(i);
        y(i) = ext_euclidean(M(i),m(i));
    end    
    
    x = 0;
    
    for i=1:len
        x = x + a(i)*M(i)*y(i);
    end
    
    x = mod(x,Mprd);
    
end

