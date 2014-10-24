%**************************************************************************
% Estimates the area of a region beneath a curve by applying the  
% Trapezoidal Rule.  (Can be used to estimate integrals.)
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
function [rsum] = int_trapezoidal(func,a,b,n)

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
        rsum = rsum + 2*func(xi);
    end
    
    rsum = ((b-a)/(2*n))*(rsum + func(a) + func(b));

end

