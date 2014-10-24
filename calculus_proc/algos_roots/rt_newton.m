%**************************************************************************
% Finds the root of a function using the Newton method, an iterative
% algorithm that uses values of f(x) to find an interval such that f(x) has
% opposite signs at the endpoints.  An initial x1 is chosen to be one
% endpoint of this interval and then the following formula is used: 
% x(n+1)=x(n) - f(x(n))/f'(x(n)), n=1,2,3
%
% @param   func   A handle to the function to find a root of.
% @param   funcp  A handle to the derivative of the above function.
% @param   a      The lower bound on the interval to find the root for.
% @param   b      The upper bound on the interval to find the root for.
% @param   tol    The (optional) tolerance: if |x(n+1)-x(n)|<tol then 
%                 x(n+1) is a root.
%
% @return  A root of the function within the interval [a,b] such that
%         |func(c)|<tol.
%
% @throw   Error if both of func(a)<0<func(b) and func(b)<0<func(a) are not 
%          true.
%
% [1] Hunt, Richard A., Calculus of a Single Variable, Harper & Row
%     1998.
%**************************************************************************
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
function [c] = rt_newton(func,funcp,a,b,tol)

    if ~exist('tol','var')
        tol = 0.0000001;
    end
    
    % ensure that a<=b
    if b<a
        t = a;
        a = b;
        b = t;
    end
    
    % ensure that one of func(a)<0<func(b) or func(b)<0<func(a)
    if ~((func(a)<0 && func(b)>0) || (func(a)>0 && func(b)<0))
       error('rt_newton: interval assumptions are not satisfied');
    end
    
    % use Newton's until a suitable solution has been found or the maximum 
    % number of iterations has been reached
    MAX_ITERS = 1000;    
    
    xn = a;

    i = 0;
    while true

        i = i+1;

        xnl = xn - func(xn)/funcp(xn);

        % has xn1 changed much?
        if abs(xnl-xn)<tol
            c = xnl;
            break;
        else
            xn = xnl;   
        end

        % if the maximum number of iterations has been reached then
        % break out
        if i>MAX_ITERS
            break;
        end

    end
    
    
    
    

    