%**************************************************************************
% Generates a exponential discrete time signal.
%
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
function [y] = exponential(a,t)

    if ~isvector(t)
        error('exponential: input t is not a vector');
    end

    y = a.^t;
end
