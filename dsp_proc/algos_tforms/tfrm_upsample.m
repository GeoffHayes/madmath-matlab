%**************************************************************************
% Performs a transformation up-sample on the independent time variable N
% and evaluates the up-sampled time data given the function func.
%
% @param   N      The time period to shift
% @param   m      The up-sampling factor.
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
function[y] = tfrm_upsample(N,m,func)
    y = zeros(length(N),1);
    
    for i=1:length(N)
        if mod(N(i),m)==0
            y(i) = func(N(i)/m);
        end
    end
end