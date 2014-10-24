%**************************************************************************
% Estimates the area of a region beneath a curve by applying the Riemann 
% Sums technique.  (Can be used to estimate integrals.)
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
function [rsum] = int_riemannsum(func,a,b,n)

    % ensure that a<=b
    if b<a
        t = a;
        a = b;
        b = t;
    end
    
    n = ceil(abs(n));
    
    delta = (b-a)/n;
    
    rsum = 0;
    xi   = a;
    
    if func(xi)==0
        xi = a + delta;
    end

    for i=1:n       
        rsum = rsum + delta*func(xi);
        xi   = xi + delta;
    end

end

