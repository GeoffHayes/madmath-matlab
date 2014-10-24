%**************************************************************************
% Miller-Rabin primality test that indicates if an integer is probable prime or
% not.
%
% @param   p   A positive integer to check for primality.
% @param   k   A postive integer indicating the number of iterations of the
%              algorithm.
%
% @return  Indicator if the integer is a probable prime (1) or not (0).
%
% [1] Stinson, Douglas R., Cryptography Theory and Practice, CRC Press
%     1995.
% [2] http://en.wikipedia.org/wiki/Miller?Rabin_primality_test
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
function [bool] = prmlty_mllrrabin(p,k)

    if ~exist('k','var')
        k = 50;
    end
    
    k = abs(round(k));
    p = cast(p,'uint64');
    
    bool = 1;
    
    if p==1
        bool = 0;
    elseif p==2 || p==3 || p==5
        bool = 1;
    elseif mod(p,2)==0
        bool = 0;
    else
        
        % write p as p-1=2^z * m where m is odd
        z    = 0;
        m    = p-1;
        mask = 1;
        
        while bitand(m,mask)~=1
            z = z + 1;
            m = bitshift(m,-1);
        end

        % iterate k times, randomly choosing an a to validate the above theorem
        a = cast(1 + floor(double(p-1)*rand(k,1)),'uint64');
        
        for i=1:k
            b = safe_moduloexp(a(i),m,p);
            
            if b==1 || b==(p-1) % equivalent to b=1modp or b=-1modp
                continue;
            end

            bool = 0;
            
            for j=0:z-1
                b = safe_moduloexp(b,2,p);
                if b==1
                    break;
                elseif b==(p-1)
                    bool = 1;
                    break;
                end
            end
            
            if ~bool
                break;
            end
        end
    end
end

