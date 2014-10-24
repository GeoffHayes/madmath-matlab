%**************************************************************************
% Sinc (sine cardinal) function defined by:
% 
%               sinc(x) = 1,              if x==0
%                       = sin(pix)/(pix), if x~=0
%
% Note that a property of the sinc function is that the set of local
% extrema of sinc(x) corresponds to its intersections with the cosine
% function.
%
% @param   x   The input array.
%
% @return  The sinc of the input.
%
% @throw   Error if the input sequence is not a vector.
%
% http://www.mathworks.com/help/signal/ref/sinc.html
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
function [y] = sinc(x)

    y       = sin(pi*x)./(pi*x);
    y(x==0) = 1;

end

