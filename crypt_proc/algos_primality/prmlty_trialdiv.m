%**************************************************************************
% Trial division primality test that indicates if an integer is a 
% prime or not.
%
% @param   p   A positive integer to check for primality.
%
% @return  Indicator if the integer is a probable prime (1) or not (0).
%
% [1] Stinson, Douglas R., Cryptography Theory and Practice, CRC Press
%     1995.
% [2] http://en.wikipedia.org/wiki/Trial_division
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
function [bool] = prmlty_trialdiv(p)

    p = cast(p,'uint64');
    
    if p==1
        bool = 0;
    elseif p==2 || p==3 || p==5
        bool = 1;
    elseif bitand(p,1)==0
        % input number is even so can't be prime
        bool = 0;
    else
        ub = cast(sqrt(double(p)),'uint64');
        fctrs = [2 3:2:ub];
        
        % divide p by all possible factors up to the square root of p
        bool = all(rem(p,fctrs));
        
    end

end

