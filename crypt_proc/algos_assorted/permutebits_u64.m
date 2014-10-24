%**************************************************************************
% Permutes the bits of a 64-bit unsigned integer according to the
% permutation matrix P.
%
% @param   x     The integer to permute.
% @param   P     The permutation matrix to permute the bits.
% @param   msb   The most significant bit (1-64) within the input integer.
%
% @return  The bit permuted integer.
%
% @throw   Error if the input integer is not a 64-bit unsigned integer.
% @throw   Error if the permutation matrix contains elements outside of the
%          range 1,2,...,64.
%
% @warning x is assumed to be written as the concatenation of 64 bits
%          numbered as follows: b1,b2,....,b64.
% [1] Stinson, Douglas R., Cryptography Theory and Practice, CRC Press
%     1995.
% [2] http://en.wikipedia.org/wiki/DES_supplementary_material
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
function [y] = permutebits_u64(x,P,msb)

    if ~exist('msb','var')
        msb = 64;
    end
    
    numelems = numel(P);
    y        = cast(0,class(x));
    
    for i=1:numelems
        y = bitor(bitshift(y,1),bitget(x,msb-P(i)+1));
    end

end