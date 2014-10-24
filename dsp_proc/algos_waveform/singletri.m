%**************************************************************************
% Generates a single triangle with peak a centred at point c  
% (i.e. f(c) = a) with slope m.
%
% @param   a   Triangle peak (height).
% @param   c   Triangle centre.
% @param   m   Triangle slope (+m on left side, -m on right side).
% @param   t   Monotonically increasing time vector to be used as input to
%              the exponential function.
%
% @throw   Error if t is not a vector.
%
% http://www.mathworks.com/help/signal/ref/square.html
% http://en.wikipedia.org/wiki/Square_wave
%**************************************************************************
% Code is written by author based upon noted references (if given), using 
% The MathWorks MATLAB function signature (if applicable) for convenience 
% only.
%
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
function [y] = singletri(a,c,m,t)

    if ~isvector(t)
        error('singletri: input t is not a vector');
    end

    y = zeros(size(t));
    
    for i=1:length(t)
       
        if t(i)<=c
            y(i) = max(0,m*t(i)-m*c+a);
        else
            y(i) = max(0,-m*t(i)+m*c+a);
        end
    end
end
