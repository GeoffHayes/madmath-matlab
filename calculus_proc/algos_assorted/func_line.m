%**************************************************************************
% Applies the linear function fx=y=mx+b to the input vector x, returning
% the resulting y vector.
%
% @param   m   The slope for the line.
% @param   b   The y-intercept.
% @param   x   The x data to evaluate y for.
%
% @return  The result of evaluating y=mx+b and the vector x (both as column
%          vectors) within a matrix as [x y].
%
% @note    If m is Inf or -Inf, then the line is considered to be vertical
%          and so the equation simply becomes x=b.
%
% [1] Hunt, Richard A., Calculus of a Single Variable, Harper & Row
%     1998.
%**************************************************************************
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
function [L] = func_line(m,b,x)

    if isinf(m)
        y=x;
        x=repmat(b,size(y));
    else
        y=m*x+b;
    end
    
    % ensure that column vectors are returned
    if size(x,1)==1
        x=x';
        y=y';
    end
    
    L = [x y];

end

