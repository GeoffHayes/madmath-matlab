%**************************************************************************
% Returns the volume of a right circular cone.
%
% @param   r   The radius of the base of the cone.
% @param   h   The height of the cone.
%
% @return  The volume, 1/3*pi*r^2*h.
%
% [1] Hunt, Richard A., Calculus of a Single Variable, Harper & Row
%     1998.
%**************************************************************************
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
function [v] = vol_cone(r,h)

    v = 1/3*pi*r^2*h;

end

