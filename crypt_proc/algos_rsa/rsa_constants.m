%**************************************************************************
% Defines constants used for the Rivest, Shamir and Adleman (RSA)
% algorithm.
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
function rsa_constants

    global MAX_PQ;  % maximum allowed product of the two primes p and q
    global NUM_BITS_PER_BYTE;  % number of bits in a byte
    global NUM_BYTES_PER_EBLK; % number of bytes per block that RSA encrypts
    global NUM_BYTES_PER_DBLK; % number of bytes per block that RSA decrypts
    global MIN_PRIME_NUMBER;   % the minium allowed prime for a p or q
    
    MAX_PQ             = cast(hex2dec('FFFFFFFFFFFFFFFF'),'uint64');
    NUM_BYTES_PER_EBLK = 1;  % allow 1 bytes or 8 bits since that is less
                             % than pq which is at least 9 bits long
    NUM_BYTES_PER_DBLK = 2;  % allow 2 bytes or 16 bits
    NUM_BITS_PER_BYTE  = 8;
    MIN_PRIME_NUMBER   = 257;   % allows for n (p*q) to be larger than the
                                % maximum allowed 8-bit unsigned integer

end

