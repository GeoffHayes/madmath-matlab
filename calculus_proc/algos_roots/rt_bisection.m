%**************************************************************************
% Finds the root of a function using the Bisection method: if a function func
% is continuous on the interval a<=x<=b, and either func(a)<0<func(b) or
% func(b)<0<func(a), then there is at least one number c such that a<c<b and
% func(c)=0.
%
% @param   func   A handle to the function to find a root of.
% @param   a      The lower bound on the interval to find the root for.
% @param   b      The upper bound on the interval to find the root for.
% @param   tol    The (optional) tolerance: if |func(c)|<tol then c is a 
%                 root.
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
function [c] = rt_bisection(func,a,b,tol)

    if ~exist('tol','var')
        tol = 0.0000001;
    end
    
    % ensure that a<=b
    if b<a
        t = a;
        a = b;
        b = t;
    end
    
    if abs(func(a))<tol 
        % choose a
        c = a;
    elseif abs(func(b))<tol
        % choose b
        c = b;
    else
    
        % ensure that one of func(a)<0<func(b) or func(b)<0<func(a)
        if ~((func(a)<0 && func(b)>0) || (func(a)>0 && func(b)<0))
           error('rt_bisection: interval assumptions are not satisfied');
        end

        % bisect until a suitable solution has been found or the maximum number
        % of iterations has been reached
        MAX_ITERS = 1000;

        i = 0;
        while true
            
            i = i+1;

            % calculate the midpoint between a and b
            c = (a+b)/2;
            
            % if the evaluation at c is close to zero then a root has been
            % found
            if abs(func(c))<tol
                break;
            else
                if (func(a)<0 && func(c)>0) || (func(a)>0 && func(c)<0)
                    % replace b with c
                    b = c;
                else
                    % replace a with c
                    a = c;
                end    
            end
            
            % if the maximum number of iterations has been reached then
            % break out
            if i>MAX_ITERS
                break;
            end

        end
    end
end

