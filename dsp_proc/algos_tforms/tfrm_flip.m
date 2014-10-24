%**************************************************************************
% Performs a transformation flip (or reversal) on the independent time 
% variable N and evaluates the flipped time data given the function func.
%
% @param   N      The time period to flip.
% @param   func   The evaluating function.
%
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
function[y] = tfrm_flip(N,func)
    Np = -N;
    y = func(Np);
end