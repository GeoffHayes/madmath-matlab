%**************************************************************************
% The p-1 factoring algorithm to determine if the input n has a factor.  If
% no factor is found, then that does not guarantee that there are none for
% the integer (just that the upper bound is not great enough).
%
% @param   n    The integer to find a factor of.
% @param   ub   The upper bound on the "search space".
%
% @return  A non-zero factor (if one exists) or zero (failure of
%          algorithm).
%
% [1] Stinson, Douglas R., Cryptography Theory and Practice, CRC Press
%     1995.
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
function [fctr] = fctr_plessone(n,ub)

    n  = uint64(round(abs(n(1))));
    ub = round(abs(ub(1)));
    
    if n<4
        fctr=n;
    elseif bitand(n,uint64(1))==0
        % n is even so the obvious factor is 2
        fctr=2;
    else
        % n is odd so apply the algorithm
        a=uint64(2);
        ap=a;
        for j=2:ub
            a = safe_moduloexp(a,j,n);
            if a==1
                % use previous a since not one
                a=ap;
                break;
            else
                % save current a to previous
                ap=a;
            end
        end
        
        fctr = euclidean(a-1,n);
        if fctr>1 && fctr<n
            % fctr is a factor of n
        else
            % fctr is not a factor
            fctr=0;
        end 
    end
end

