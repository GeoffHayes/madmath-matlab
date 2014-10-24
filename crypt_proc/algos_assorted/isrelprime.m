%**************************************************************************
% Determines whether all pairs of integers are relatively prime. i.e. the 
% GCD of any two integers is one.
%
% @param   x   Vector of integers to determine if all are pairwise
%              relatively prime.
%
% @return  Flag indicating if all integers are relatively prime (1) or
%          not (0).
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
function [bool] = isrelprime(x)

    if isscalar(x)
        % x is relatively prime to itself
        bool = 1;
    else
        % x is a vector or matrix so compare all pairs
        bool     = 1;
        numelems = numel(x);
        for i=1:numelems
            for j=i+1:numelems
                if euclidean(x(i),x(j))~=1
                    bool = 0;
                    break;
                end
            end
            if ~bool
                break;
            end
        end
    end
end

