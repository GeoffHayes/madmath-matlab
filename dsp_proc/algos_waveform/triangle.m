%**************************************************************************
% Generates a triangle waveform with period 2pi.
%
% @param   t   Monotonically increasing time vector to be used as input to
%              the triangle function.
%
% @throw   Error if t is not a vector.
%
% http://www.mathworks.com/help/signal/ref/sawtooth.html
% http://en.wikipedia.org/wiki/Triangle_wave
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
function [y] = triangle(t)

    if ~isvector(t)
        error('triangle: input t is not a vector');
    end
    
    % default the period to 2pi
    a = 2*pi;
    
    y = abs(2*(t./a - floor(0.5+t./a)));
    
    % alternative
    % y = 2*(t./a - floor(0.5+t./a))*(-1).^floor(t./a+0.5)
end

% 50 Hz triangle displaying only 1/5 of a second (see T)
% T = 10*(1/50);
% Fs = 1000;
% dt = 1/Fs;
% t = 0:dt:T-dt;
% x = triangle(2*pi*50*t);
% plot(t,x)