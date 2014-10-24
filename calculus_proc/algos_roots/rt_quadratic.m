%**************************************************************************
% Finds the roots of a quadratic polynomial.
%
% @param   a   The scalar for the second-degree term of the polynomial.
% @param   b   The scalar for the first-degree term of the polynomial.
% @param   c   The constant term of the polynomial.
%
% @return  A 2x1 vector of real or complex roots.
%
% [1] Hunt, Richard A., Calculus of a Single Variable, Harper & Row
%     1998.
%**************************************************************************
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
function [roots] = rt_quadratic(a,b,c)

    roots = zeros(2,1);
    
    t        = sqrt(b^2 - 4*a*c);
    roots(1) = (-b + t)/(2*a);
    roots(2) = (-b - t)/(2*a);  

end

