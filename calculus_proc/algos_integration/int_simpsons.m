%**************************************************************************
% Estimates the area of a region beneath a curve by applying Simpson's  
% Rule.  (Can be used to estimate integrals.)
%
% @param   func   A handle to the function of a curve.
% @param   a      The lower bound on the interval of region of interest.
% @param   b      The upper bound on the interval of region of interest.
% @param   n      The number of sub-intervals within [a,b].
%
% @return  The estimated area.
%
% [1] Hunt, Richard A., Calculus of a Single Variable, Harper & Row
%     1998.
%**************************************************************************
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
function [rsum] = int_simpsons(func,a,b,n)

    % ensure that a<=b
    if b<a
        t = a;
        a = b;
        b = t;
    end
    
    n = ceil(abs(n));
    
    % note that x0==a and xn==b
    delta = (b-a)/(n);
    
    rsum = 0;
    xi   = a;
    
    for i=1:n-1   
        xi   = xi + delta;
        if mod(i,2)==0
            c = 2;
        else
            c = 4;
        end
        rsum = rsum + c*func(xi);
    end
    
    rsum = ((b-a)/(3*n))*(rsum + func(a) + func(b));

end

