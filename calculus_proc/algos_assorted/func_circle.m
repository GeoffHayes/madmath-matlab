%**************************************************************************
% Applies the equation of the circle (x-h)^2 + (y-k)^2 = r^2, returning the
% resulting x and y vectors.
%
% @param   h   The x coordinate for the centre of the circle.
% @param   k   The y coordinate for the centre of the circle.
% @param   r   The radius of the circle.
% @param   d   The optional resolution along the x-axis (0.1 is the
%              default).
%
% @return  The result of evaluating the above equation and the vector x 
%          within a matrix as [x y1 y2] where y1 is a vector of points for 
%          the upper half of the circle, and y2 is a vector of points for 
%          the lower half of the circle.
%
% [1] Hunt, Richard A., Calculus of a Single Variable, Harper & Row
%     1998.
%**************************************************************************
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
function [C] = func_circle(h,k,r,d)

    if ~exist('d','var')
        d=0.1;
    end

    r = abs(r);
    
    % create the vector of points along the x-axis
    x = (-r:d:r)';
    
    C = zeros(length(x),3);
    
    % update C according to the equation of a circle
    C(:,1) = x;
    C(:,2) = sqrt(r^2 - x.^2);
    C(:,3) = -C(:,2);
    
    % apply the shifts
    C(:,1) = C(:,1)+h;
    C(:,2) = C(:,2)+k;
    C(:,3) = C(:,3)+k;
end

