%**************************************************************************
% Performs a safe modulo operation against a positive integer raised to
% some positive exponent: mod(a^k,m)
%
% @param   a   A positive integer to perform the modulo operation against.
% @param   k   A positive integer exponent for a.
% @param   m   A postive integer to reduce a^k by modulo m.
%
% @return  The result of mod(a^k,m)
%
% @throw   Error if integer overflow occurs due to the multiplication.
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
function [r] = safe_moduloexp(a,k,m)

    a = cast(abs(round(a)),'uint64');
    k = cast(abs(round(k)),'uint64');
    m = cast(abs(round(m)),'uint64');
    
    if m==1
        r = 0;
    else
        % reduce a modulo m
        a = mod(a,m);
    
        if k==0
            r = mod(1,m);
        elseif k==1
            r = a;
        elseif k==2
            % will a^2 cause an overflow?
            r=validate_moduloprod(a,a,m);
        elseif mod(k,2)==0
            r = safe_moduloexp(a,k/2,m);
            % will ap^2 cause an overflow?
            r=validate_moduloprod(r,r,m);
        else
            r = safe_moduloexp(a,k-1,m);
            % will a*ap cause an overflow?
            r=validate_moduloprod(r,a,m);
        end
    end
end

%**************************************************************************
% Validates whether the product of two integers can proceed i.e. doesn't
% cause an integer overflow exception.  If the exception will occur then
% attempts are made to mitigate the overflow and return the result.
%
% @param   a   A positive integer to perform the modulo operation against.
% @param   b   A positive integer exponent for a.
% @param   m   A postive integer to reduce a^k by modulo m.
%
% @return  The result of mod(a*b,m)
%
% @throw   Error if integer overflow occurs due to the multiplication.
%**************************************************************************
% Code is written by author based upon noted references, using the
% MATLAB signature (if applicable) for convenience only.  This software 
% should in no way be considered that as a version (official or otherwise) 
% provided by Mathworks.
%
% Geoff Hayes 2014
%**************************************************************************
function [r] = validate_moduloprod(a,b,m)

    MAX_UINT64 = cast(hex2dec('FFFFFFFFFFFFFFFF'),'uint64');

    if MAX_UINT64/b < a
        
        % cast as 128-bit unsigned ints
        a128=uint128(a);
        b128=uint128(b);
        m128=uint128(m);

        % perform the multiplication
        b128=a128*b128;
        
        % determine the msb for each integer and compute the maximum number
        % of bitshifts for m128
        mmsb=m128.msb;
        bmsb=b128.msb;
        k   =abs(mmsb-bmsb)+1;
        
        m128p=m128.shift(k);
        
        % keep subtracting out large multiples of m (m128p) until the upper
        % 64 bits of b128 are zero (or until k cannot be reduced further)
        while k>0 && b128.msb64~=0
           if m128p>b128
               k=k-1;
               m128p=m128p.shift(-1);
           end
           while b128>=m128p
               b128=b128-m128p;
           end
        end
        
        % if the upper 64 bits of b128 are zero, then use usual mod to
        % reduce the lower 64 bits of b128 modulo m
        if b128.msb64==0
            r=mod(b128.lsb64,m);
        else
            error('validate_moduloprod: cannot reduce product modulo %d',m);
        end
    else
        r = mod(a*b,m);
    end
end

