%**************************************************************************
% Applies the equation of the ellipse ((x-h)^2)/a^2 + ((y-k)^2)/b^2 = 1, 
% returning the resulting x and y vectors.
%
% @param   h   The x coordinate for the centre of the ellipse.
% @param   k   The y coordinate for the centre of the ellipse.
% @param   a   The distance from the centre of ellipse to one of (the two)
%              points of interesection of the ellipse with the x-axis.
% @param   a   The distance from the centre of ellipse to one of (the two)
%              points of interesection of the ellipse with the y-axis.
% @param   d   The optional resolution along the x-axis (0.1 is the
%              default).
%
% @return  The result of evaluating the above equation and the vector x 
%          within a matrix as [x y1 y2] where y1 is a vector of points for 
%          the upper half of the ellipse, and y2 is a vector of points for 
%          the lower half of the ellipse.
%
% [1] Hunt, Richard A., Calculus of a Single Variable, Harper & Row
%     1998.
% [2] http://en.wikipedia.org/wiki/Ellipse
%
%**************************************************************************
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
function [C] = func_ellipse(h,k,a,b,d)

    if ~exist('d','var')
        d=0.1;
    end

    a = abs(a);
    b = abs(b);
    
    % create the vector of points along the x-axis, shifted by h
    x = (-a:d:a)';
    
    C = zeros(length(x),3);
    
    % update C according to the equation of a circle
    C(:,1) = x;
    C(:,2) = sqrt(1 - (x.^2)./(a^2))*b;
    C(:,3) = -C(:,2);
    
    % apply the shifts
    C(:,1) = C(:,1)+h;
    C(:,2) = C(:,2)+k;
    C(:,3) = C(:,3)+k;
end

