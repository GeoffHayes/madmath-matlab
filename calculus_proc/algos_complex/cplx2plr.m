%**************************************************************************
% Converts a complex number to its equivalent in polar coordinates.
%
% @param   cplx   The complex number
% @return  The polar coordinate.
%
% [1] Stewart, James, Multivariable Calculus, Brooks/Cole Publishing
%     Company, 1991.
%**************************************************************************
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
function [r,theta] = cplx2plr(c)

    a = real(c);
    b = imag(c);
    
    r = sqrt(a^2 + b^2);
    theta = atan2(b,a);

end