%**************************************************************************
% Demoivre's theorem - if n is a positive integer and z=r(cosa + isina)
% then z^n = r^n(cosna + isinna).
%
% @param   cplx   The complex number.
% @param   n      The exponent n.
%
% @return  The comple number to the power of n.
%
% [1] Stewart, James, Multivariable Calculus, Brooks/Cole Publishing
%     Company, 1991.
%**************************************************************************
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
function [cplx] = thrm_demoivre(cplx,n)

    n = abs(ceil(n));

    % convert to polar coordinates
    [r,theta] = cplx2plr(cplx);
    
    cplx      = r^n*(cos(n*theta) + 1i*sin(n*theta));

end