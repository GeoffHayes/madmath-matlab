%**************************************************************************
% Performs a transformation down-sample on the independent time variable N
% and evaluates the down-sampled time data given the function func.
%
% @param   N      The time period to shift
% @param   m      The downsampling factor.
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
function[y] = tfrm_downsample(N,m,func)
    Np = m*N;
    y = func(Np);
end