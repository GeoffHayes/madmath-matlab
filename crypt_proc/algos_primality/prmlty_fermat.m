%**************************************************************************
% Fermat primality test that indicates if an integer is probable prime or
% not.
%
% @param   p   A positive integer to check for primality.
% @param   k   A postive integer indicating the number of iterations of the
%              algorithm.
%
% @return  Indicator if the integer is a probable prime (1) or not (0).
%
% [1] http://en.wikipedia.org/wiki/Fermat_primality_test
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
function [bool] = prmlty_fermat(p,k)

    if ~exist('k','var')
        k = 50;
    end
    
    k = abs(round(k));
    p = cast(p,'uint64');
    
    % Fermat's little theorem states that if p is prime and 1<=a<p then,
    % a^(p-1)=1 mod p
    
    bool = 1;
    
    if p==1
        bool = 0;
    elseif p==2 || p==3 || p==5
        bool = 1;
    elseif mod(p,2)==0
        bool = 0;
    else
        % iterate k times, randomly choosing an a to validate the above theorem
        a = cast(1 + floor(double(p-1)*rand(k,1)),'uint64');
        
        for i=1:k
            if ~(safe_moduloexp(a(i),p-1,p)==1)
                bool = 0;
                break;
            end
        end
    end
end

