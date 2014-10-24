%**************************************************************************
% Executes the extended Euclidean algorithm to compute the inverse of
% integer b modulo m.
%
% @param   b   Integer to find the inverse for.
% @param   m   The size of the symbol space.
%
% @return  The inverse (if it exists) of b.
%
% @throw   Error if b has no inverse modulo m.
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
function [ib] = ext_euclidean(b,m)

    ctype = class(b);
    
    b = double(abs(round(b(1))));
    m = double(abs(round(m(1))));
    
    if b>=m
        b = mod(b,m);
    end
    
    m0 = m;
    b0 = b;
    t0 = 0;
    t  = 1;
    q  = floor(m0/b0);
    r  = m0 - q*b0;
      
    while r

        temp = t0-q*t;
        if temp>=0
            temp = mod(temp,m);
        elseif temp<0
            temp = m-mod(-temp,m);
        end
        
        t0 = t;
        t  = temp;
        m0 = b0;
        b0 = r;
        q = floor(m0/b0);
        r = m0 - q*b0;

    end
    
    if b0~=1
        error('ext_euclidean: b (%d) has no inverse modulo m (%d)',b,m);
    else
        ib = mod(t,m);
    end
    
    ib = cast(ib,ctype);
        
end

