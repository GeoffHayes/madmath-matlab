%**************************************************************************
% Returns the volume of a right circular cylinder.
%
% @param   r   The radius of the base (or top) of cylinder.
% @param   h   The height of the cylinder.
%
% @return  The volume, pi*r^2*h.
%
% [1] Hunt, Richard A., Calculus of a Single Variable, Harper & Row
%     1998.
%**************************************************************************
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
function [v] = vol_cylinder(r,h)

    v = pi*r^2*h;

end

