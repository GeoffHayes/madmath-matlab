%**************************************************************************
% Class definition for an 128-bit unsigned integer.
%
%**************************************************************************
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
classdef uint128
    % Class to represent a 128-bit unsigned integer.
    
    properties
        % Most significant (upper) 64-bits (MSBs) of the unsigned integer.
        msb64;        
        % Least significant (lower) 64-bits (LSBs) of the unsigned integer.
        lsb64;
    end
    
    properties (Constant, Access=private)
       % mask of all zeros except for bit 64
       MASK_BIT64_UINT64 = cast(hex2dec('8000000000000000'),'uint64');
       % mask of all zeros except for bit 1
       MASK_BIT1_UINT64  = uint64(1);
       % number of bits in a 64-bit unsigned integer
       NUM_BITS_UINT64   = 64;
       % maximum allowed 64-bit unsigned integer
       MAX_UINT64        = cast(hex2dec('FFFFFFFFFFFFFFFF'),'uint64');
       % number of bits in a 128-bit unsigned integer
       NUM_BITS_UINT128  = 128;
       % mask of all zeros
       MASK_0_UINT64     = 0;
    end  

    methods
       
        % 128-bit unsigned integer class constructor.
        %
        % @param   y   The optional input parameter to initialize the
        %              128-bit unsigned integer.
        %
        % @return  The new 128-bit unsigned integer.
        %
        function [x] = uint128(varargin)
            x.msb64 = uint128.MASK_0_UINT64;
            x.lsb64 = x.msb64;
            
            if nargin==1
            
                y     = varargin{1};
                ctype = class(y);

                switch ctype

                    case 'uint128'
                        x.msb64 = y.msb64;
                        x.lsb64 = y.lsb64;
                    case {'uint64','uint32','uint16','uint8'}
                        x.lsb64 = cast(y,'uint64');
                    case {'int64','int32','int16','int8'}
                        x.lsb64 = cast(y,'uint64');
                    case {'single','double'}
                        x.lsb64 = cast(y,'uint64');
                    otherwise
                        error('uint128: invalid input to the constructor');
                end
                
            end

        end 
        
        % 128-bit unsigned integer addition operator.
        %
        % @param  x   The uint128 on the left side of the addition.
        % @param  y   The uint128 on the right side of the addition.
        %
        % @return  The addition of x and y (x+y).
        %
        % @throw  Error if the inputs are not of type uint128.
        %
        % @warning If x+y>=2^128, then x+y=2^128-1.
        function [r] = plus(x,y)

            r   = uint128;
            
            % add the LSB bits then the MSB bits
            cry = false;
            xuint64 = x.lsb64;
            yuint64 = y.lsb64;
            ruint64 = r.lsb64;
            for j=1:2

                % reset the mask
                msk = uint128.MASK_BIT1_UINT64;                
                
                % if one of x or y is zero, then just copy the other and
                % skip to the next j iteration
                if (xuint64==0 || yuint64==0)
                    if xuint64==0
                        ruint64=yuint64;
                    else
                        ruint64=xuint64;
                    end
                    
                    if j==1
                        r.lsb64 = ruint64;

                        xuint64 = x.msb64;
                        yuint64 = y.msb64;
                        ruint64 = r.msb64;
                    else
                        r.msb64 = ruint64;
                        
                        if cry
                            r.msb64 = r.msb64+1;
                            cry = 0; 
                        end
                    end
                    
                    continue;
                end
                
                for i=1:uint128.NUM_BITS_UINT64

                    % apply the carry-over
                    if cry
                        ruint64 = bitor(ruint64,msk);
                        cry     = false;
                    end

                    if bitand(xuint64,msk) 
                        if bitand(yuint64,msk)
                            % 1 + 1 = 0 so do nothing
                            % set the carry-over flag to true
                            cry = true;
                        else
                            if bitand(ruint64,msk)
                                % 1 + 1 = 0
                                ruint64 = bitxor(ruint64,msk);
                                cry     = true;
                            else
                                % 1 + 0 = 1
                                ruint64 = bitor(ruint64,msk);
                            end
                        end
                    elseif bitand(yuint64,msk)
                        if bitand(ruint64,msk)
                            % 1 + 1 = 0
                            ruint64 = bitxor(ruint64,msk);
                            cry     = true;
                        else
                            % 1 + 0 = 0
                            ruint64 = bitor(ruint64,msk);
                        end  
                    end

                    msk = bitshift(msk,1);
                end
                
                
                if j==1
                    r.lsb64 = ruint64;
                    
                    xuint64 = x.msb64;
                    yuint64 = y.msb64;
                    ruint64 = r.msb64;
                else
                    r.msb64 = ruint64;
                end
            end
            
            % if a final bit is to be carried over, then the maxium has
            % been reached so force it
            if cry
               r.lsb64 = uint128.MAX_UINT64;
               r.msb64 = r.lsb64;
            end
        end 
        
        % 128-bit unsigned integer subtraction operator.
        %
        % @param  x   The uint128 on the left side of the subtraction.
        % @param  y   The uint128 on the right side of the subtraction.
        %
        % @return The difference between x and y (x-y).
        %
        % @throw  Error if the inputs are not of type uint128.
        %
        % @warning If x-y<0, then x-y=0.        
        function [r] = minus(x,y)

            r = uint128;
            
            % can the subtraction proceed?
            if x.msb64<y.msb64 || (x.msb64==y.msb64 && x.lsb64<y.lsb64)
                % y is greater than x so keep the result as zero
            elseif x.msb64==y.msb64 && x.lsb64==y.lsb64
                % x and y are identical so result is zero
            elseif x.msb64==0 && y.msb64==0
                % upper bits are all zeros for both x and y, hence need
                % only subtract the lower bits
                r.lsb64 = x.lsb64-y.lsb64;
            else
            
                % subtract the LSB bits then the MSB bits
                brw     = false;
                xuint64 = x.lsb64;
                yuint64 = y.lsb64;
                ruint64 = r.lsb64;
                for j=1:2
                    % reset the mask
                    msk = cast(1,'uint64');

                    for i=1:uint128.NUM_BITS_UINT64

                        % apply the borrow if possible
                        if brw 
                            if bitand(xuint64,msk)
                                xuint64 = bitxor(xuint64,msk);
                                brw     = false;
                            else
                                % borrow again 
                                xuint64 = bitor(xuint64,msk);
                            end
                        end

                        if bitand(xuint64,msk) 
                            if bitand(yuint64,msk)
                                % do nothing as 1 - 1 = 0
                            else
                                % 1 - 0 = 1
                                ruint64 = bitor(ruint64,msk);
                            end
                        elseif bitand(yuint64,msk)
                            % 0 - 1 = 1 
                            ruint64 = bitor(ruint64,msk);
                            brw     = true;
                        end
                        
                        % do nothing for 0 - 0 = 0
                        msk = bitshift(msk,1);
                    end


                    if j==1
                        r.lsb64 = ruint64;

                        xuint64 = x.msb64;
                        yuint64 = y.msb64;
                        ruint64 = r.msb64;
                    else
                        r.msb64 = ruint64;
                    end
                end   
            end
        end  
        
        % 128-bit unsigned less than comparison operator.
        %
        % @param  x   The uint128 on the left side of the operator.
        % @param  y   The uint128 on the right side of the operator.
        %
        % @return True if x<y and false otherwise.
        %
        % @throw  Error if the inputs are not of type uint128.
        %
        function [bool] = lt(x,y)

            bool =  y.msb64>x.msb64 || ...
                (y.msb64==x.msb64 && y.lsb64>x.lsb64);
        end
            
        % 128-bit unsigned greater than comparison operator.
        %
        % @param  x   The uint128 on the left side of the operator.
        % @param  y   The uint128 on the right side of the operator.
        %
        % @return True if x>y and false otherwise.
        %
        % @throw  Error if the inputs are not of type uint128.
        %
        function [bool] = gt(x,y)  

            bool =  y.msb64<x.msb64 || ...
                (y.msb64==x.msb64 && y.lsb64<x.lsb64);
        end
        
        % 128-bit unsigned less than or equal to comparison operator.
        %
        % @param  x   The uint128 on the left side of the operator.
        % @param  y   The uint128 on the right side of the operator.
        %
        % @return True if x<=y and false otherwise.
        %
        % @throw  Error if the inputs are not of type uint128.
        %
        function [bool] = le(x,y)            
            
            bool = ~gt(x,y);
        end
            
        % 128-bit unsigned greater than comparison operator..
        %
        % @param  x   The uint128 on the left side of the operator.
        % @param  y   The uint128 on the right side of the operator.
        %
        % @return True if x>=y and false otherwise.
        %
        % @throw  Error if the inputs are not of type uint128.
        %
        function [bool] = ge(x,y)  
            
            bool = ~lt(x,y);
        end    
        
        % 128-bit unsigned equality operator.
        %
        % @param  x   The uint128 on the left side of the operator.
        % @param  y   The uint128 on the right side of the operator.
        %
        % @return True if x>=y and false otherwise.
        %
        % @throw  Error if the inputs are not of type uint128.
        %
        function [bool] = eq(x,y)              
            
            bool = x.msb64==y.msb64 && x.lsb64==y.lsb64;
        end  
        
        % 128-bit unsigned not equal operator.
        %
        % @param  x   The uint128 on the left side of the operator.
        % @param  y   The uint128 on the right side of the operator.
        %
        % @return True if x>=y and false otherwise.
        %
        % @throw  Error if the inputs are not of type uint128.
        %
        function [bool] = ne(x,y)               
            
            bool = ~eq(x,y);
        end          

        % 128-bit unsigned multiplication operator.
        %
        % @param  x   The uint128 on the left side of the operator 
        %             (multiplicand).
        % @param  y   The uint128 on the right side of the operator
        %             (multiplier).
        %
        % @return The product x*y.
        %
        % @warning If x*y>=2^128, then x*y=2^128-1.
        %
        % @throw  Error if the inputs are not of type uint128.
        %
        function [r] = mtimes(x,y)
      
            r   = uint128;
            
            % is the multiplication necessary?
            if (x.msb64==0 && x.lsb64==0) || (y.msb64==0 && y.lsb64==0)
                % one of x or y is zero, so the product is zero
            else
                % temp variable
                tmp = uint128;
                
                xluint64 = x.lsb64;
                yluint64 = y.lsb64;
                
                xmuint64 = x.msb64;
                ymuint64 = y.msb64;
                
                % initialize temp
                tmp.lsb64 = xluint64;
                tmp.msb64 = xmuint64;
                
                msk   = uint128.MASK_BIT1_UINT64;
                msk1  = uint128.MASK_BIT1_UINT64;
                msk64 = uint128.MASK_BIT64_UINT64;

                if bitand(yluint64,msk)
                    % the first bit of y is a one so copy x into the result
                    r.lsb64 = xluint64;
                    r.msb64 = xmuint64;
                    
                    % remove the bit
                    yluint64 = bitxor(yluint64,msk);
                end
                
                msk = bitshift(msk,1);
                cry = 0;
                
                % do the multiplication
                i=1;
                while ~(yluint64==0 && ymuint64==0)
                    
                    i = i+1;
                    
                    if i==(uint128.NUM_BITS_UINT64+1)
                        % reset the mask for the msb
                        msk = uint128.MASK_BIT1_UINT64;
                    end
                    
                    % shift the temp variables
                    if bitand(tmp.lsb64,msk64)
                        % the bitshift will result in a carry over to the
                        % upper set of bits
                        cry = 1;
                    end
                    
                    tmp.lsb64 = bitshift(tmp.lsb64,1);
                    
                    if bitand(tmp.msb64,msk64)
                        % the shift will result in an overflow so max the 
                        % 128-bit integers and break out
                        r.lsb64 = uint128.MAX_UINT64;
                        r.msb64 = uint128.MAX_UINT64;
                        break;
                    end
                    
                    tmp.msb64 = bitshift(tmp.msb64 ,1);
                    
                    if cry
                        % carry over so set the first bit to one
                        tmp.msb64  = bitxor(tmp.msb64 ,msk1);
                        cry        = 0;
                    end
                    
                    if i<=uint128.NUM_BITS_UINT64 && bitand(yluint64,msk)
                        % sum r with the temp
                        r = r+tmp;
                        
                        % remove the bit
                        yluint64 = bitxor(yluint64,msk);
 
                    elseif i>uint128.NUM_BITS_UINT64 && bitand(ymuint64,msk)
                        % sum r with the temp
                        r = r+tmp;    
                        
                        % remove the bit
                        ymuint64 = bitxor(ymuint64,msk);                      
                    end
                    
                    msk = bitshift(msk,1);
                end
            end
        end
        
        % Returns the most-significant bit position within x i.e. the
        % left-most one of the 128-bit unsigned integer, where x is bit
        % aligned as: b128,b127,b126,....,b2,b1.
        %
        % @param   x   The integer to find the most-significan bit.
        %
        % @return  The bit position (1-128) of the most significant bit.
        %
        %
        function [r] = msb(x)
            
            if ~strcmpi(class(x),'uint128')
                x = uint128(x);
            end
            
            r    = 128;
            mask = uint128.MASK_BIT64_UINT64;
            
            if x.msb64~=0
                if bitand(x.msb64,mask)==0
                    while true
                        mask = bitshift(mask,-1);
                        r    = r-1;
                        if bitand(x.msb64,mask)
                            break;
                        end
                    end
                end
            else
                r = uint128.NUM_BITS_UINT64;
                if bitand(x.lsb64,mask)==0
                    while true
                        mask = bitshift(mask,-1);
                        r    = r-1;
                        if bitand(x.lsb64,mask)
                            break;
                        end
                    end
                end
            end
        end
        
        % Executes a left (positive) or right (negative) bit shift on the
        % 128-bit unsigned integer.
        %
        % @param   x   The integer to bit shift.
        %
        % @return  The bit position (1-128) of the most significant bit.
        %
        %
        function [x] = shift(x,k)
  
            if ~strcmpi(class(x),'uint128')
                x = uint128(x);
            end
            
            k = round(k(1));
            
            if k>0
                % bit shifting to the left
                for i=1:k
                   % need to keep track of bits that are carried over from
                   % the lsb to the msb
                   cry=(bitand(x.lsb64,uint128.MASK_BIT64_UINT64)~=0);
                   x.lsb64 = bitshift(x.lsb64,1);
                   x.msb64 = bitshift(x.msb64,1);
                   if cry
                       x.msb64=bitor(x.msb64,uint128.MASK_BIT1_UINT64);
                   end
                end
            elseif k<0
                % bit shifting to the right
                k = abs(k);
                for i=1:k
                   % need to keep track of bits that are carried over from
                   % the msb to the lsb
                   cry=(bitand(x.msb64,uint128.MASK_BIT1_UINT64)~=0);
                   x.msb64 = bitshift(x.msb64,-1);
                   x.lsb64 = bitshift(x.lsb64,-1);
                   if cry
                       x.lsb64=bitor(x.lsb64,uint128.MASK_BIT64_UINT64);
                   end
                end
            end
        end
        
        % Performs modulo m arithmetic for the integer to some exponent b.
        %
        % @param   x   The integer to be raised to some exponent and the
        %              resulting product reduced modulo m.
        % @param   k   The exponent to raise x to.
        % @param   m   The modulus.
        %
        % @return  The result of raising x to the exponent k and reduced by
        %          the modulus m.
        %
        % @throw   Error if the inputs are not of type uint128.
        %
        %
        function [r] = modexp(x,k,m)
            
            mask = uint128.MASK_BIT1_UINT64;

            r = uint128;
            
            % if m==1, then there is nothing to do
            if m.lsb64~=0 || m.msb64~=0
            
                % reduce x modulo m
                if x.msb64~=0 || m.msb64~=0
                    while x>=m
                        x=x-m;
                    end
                else
                    x.lsb64=mod(x.lsb64,m.lsb64);
                end

                if k.lsb64==0 && k.msb64==0
                    % exponent is zero, so set the result to one
                    r = uint128(1);
                elseif k.lsb64==1 && k.msb64==0
                    % exponent is one, and since x has already been reduced
                    % modluo m, just assign to the output r
                    r = x;
                elseif k.lsb64==2 && k.msb64==0
                    % square x and reduce modulo m
                    r = x*x;
                    if r.msb64~=0 || m.msb64~=0
                        while r>=m
                            r=r-m;
                        end
                    else
                        r.lsb64=mod(r.lsb64,m.lsb64);
                    end
                elseif bitand(k.lsb64,mask)==0
                    % k is even so halve it and recurse
                    
                    % carry over a bit from the msb to the lsb?
                    cry = bitand(k.msb64,mask);
                    
                    % halve the integer by shifting all bits to right
                    k.msb64 = bitshift(k.msb64,-1);
                    k.lsb64 = bitshift(k.lsb64,-1);
                    if cry
                        % apply the carry over from the msb
                        umsk = uint128.MASK_BIT64_UINT64;
                        k.lsb64 = bitand(k.lsb64,umsk);
                    end
                    
                    % recurse
                    r = modexp(x,k,m);
                    
                    % square and reduce modulo m
                    r = r*r;
                    if r.msb64~=0 || m.msb64~=0
                        while r>=m
                            r=r-m;
                        end
                    else
                        r.lsb64=mod(r.lsb64,m.lsb64);
                    end
                else
                    % k is odd, so subtract one and recurse
                    r = modexp(x,k-uint128(1),m);
                    
                    % multiply by x and reduce modulo m
                    r = r*x;
                    if r.msb64~=0 || m.msb64~=0
                        while r>=m
                            r=r-m;
                        end
                    else
                        r.lsb64=mod(r.lsb64,m.lsb64);
                    end                   
                end
            end
        end
        
        % Displays the 128-bit unsigned integer to the console.
        %
        % @param   x   The integer to display/write to the console.
        %
        function display(x)
            
            uint128_constants;
            
            global DISPLAY_BINRARY;
            global POWERS2_0_127;
           
            % write the variable name first
            fprintf('\n%s = \n\t',inputname(1));
  
            if DISPLAY_BINRARY
                % write out the binary starting with the MSB
                msk = uint128.MASK_BIT64_UINT64;
                for i=1:uint128.NUM_BITS_UINT64
                   if bitand(x.msb64,msk)
                       fprintf('1');
                   else
                       fprintf('0');
                   end
                   msk = bitshift(msk,-1);
                end

                % now write out the binary LSB
                msk = uint128.MASK_BIT64_UINT64;
                for i=1:uint128.NUM_BITS_UINT64
                   if bitand(x.lsb64,msk)
                       fprintf('1');
                   else
                       fprintf('0');
                   end
                   msk = bitshift(msk,-1);
                end  
            else
                % write out the decimal value
                % extract its bits to a vector
                vctr = cast(zeros(uint128.NUM_BITS_UINT128,1),'uint16');
                
                % do lsb 
                msk = uint128.MASK_BIT1_UINT64;  
                for i=1:uint128.NUM_BITS_UINT64
                    if bitand(x.lsb64,msk)
                        vctr(i) = 1;
                    end
                    msk = bitshift(msk,1);
                end
                
                % do msb 
                msk = uint128.MASK_BIT1_UINT64;  
                for i=1:uint128.NUM_BITS_UINT64
                    if bitand(x.msb64,msk)
                        vctr(i+uint128.NUM_BITS_UINT64) = 1;
                    end
                    msk = bitshift(msk,1);
                end
                
                [c] = size(POWERS2_0_127,2);
                
                % create the output string
                str = blanks(c);
                
                % iterate over all powers of 2
                cry   = 0;
                nzidx = c;
                for i=c:-1:1
                    % apply the vector for 2 to the power of i plus the
                    % carry over
                    valstr = num2str(sum(vctr.*POWERS2_0_127(:,i)) + cry);
                    
                    cry = 0;
                    
                    % grab the last character of the string
                    str(i) = valstr(end);
                    
                    if ~strcmpi(str(i),'0')
                        nzidx=i;
                    end
                    
                    % get the carry over
                    if length(valstr)>1
                        cry = cast(str2num(valstr(1:end-1)),'uint16'); %#ok<ST2NM>
                    end
                end  
                fprintf('%s',str(nzidx:end));
            end
            fprintf('\n');
        end
    end
end

