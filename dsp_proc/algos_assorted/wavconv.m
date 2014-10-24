%**************************************************************************
% Performs a convolution on two waves.  Note that the result should be
% identical to that which is returned by MATLAB's conv function.
%
% @param   u       Input vector u to convolve with v.
% @param   v       Input vector v to convole with u.
%
% @return  The convolution of u with v according to 
%                   w[n] = SUM{k=0..p+q+1}u(k)v(n-k)
%          where p and q are the lengths of u and v respectively.
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
function [w] = wavconv(u,v)

    if ~isvector(u) || ~isvector(v)
        error('wavconv: u and v must be vectors');
    end
    
    p = length(u);
    q = length(v);
    
    w = zeros(p+q-1,1);
    
    for j=0:p+q-2 % -2 since starting at 0
        
       sum = 0;
       for k=0:p+q-2 % -2 since starting at 0
           if k<p && (j-k)>=0 && (j-k)<q
               sum = sum + u(k+1)*v(j-k+1);
           end
       end
       w(j+1) = sum;
    end
end