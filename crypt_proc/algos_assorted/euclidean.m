%**************************************************************************
% Executes the Euclidean algorithm to compute the greatest common divisor
% (GCD) of two positive integers.
%
% @param   x   First integer of pair to determine GCD.
% @param   y   Second integer of pair to determine GCD.
%
% @return  The GCD of the two integers.
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
function [gcd] = euclidean(x,y)

    x = abs(round(x(1)));
    y = abs(round(y(1)));
    
    if x==y
        gcd=x;
    else
        if x>y
            r0=x;
            r1=y;
        else
            r0=y;
            r1=x;
        end

        while r1
           
            t  = r1;
            r1 = mod(r0,r1);
            r0 = t;
            
        end
        
        gcd=r0;
    end
end

