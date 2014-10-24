%**************************************************************************
% Executes the Sieve of Eratosthenes in finding all prime numbers less than
% or equal to the input positive integer.
%
% @param   n   A positive integer upper bound on all prime integers to 
%              query for.
%
% @return  A list of prime numbers less than or equal to n.
%
% [1] http://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
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
function [lprimes] = sieve_erat(n)

    n = abs(floor(n(1)));
    
    if strcmpi(class(n),'double')~=1
        n = double(n);
    end

    % initialize the list to include the prime numbers 1,2 and 3 and all
    % odd numbers greater than 3 up to n
    lprimes = [2 3:2:n];

    % mark all non-primes (which are multiples of primes) as zero for all 
    % positive integers starting from 3
    len = length(lprimes);
    
    % need only consider those odd numbers up to the square root of n since
    % if p is prime and p is greater than the square root of n:
    %                   p > sqrt(n) ==> p^2 > n
    % and so can be excluded from our search
    ub  = sqrt(n);
    for i=3:2:ub
        
        p = lprimes((i+1)/2);
        
        % if non-zero (i.e. not marked as non-prime) then consider all
        % multiples of it from p^2 as multiples of p(p-1), p(p-2), etc.
        % would already have been accounted for by those primes less than p
        if p
            lprimes(((p^2+1)/2):p:len) = 0; 
        end
        
    end
    
    lprimes = cast(lprimes(lprimes>0),'uint64');

end

