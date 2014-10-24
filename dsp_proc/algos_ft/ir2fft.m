%**************************************************************************
% Performs the inverse radix-2 Fast Fourier Transform (FFT) on the discrete 
% sequence of frequency-domain sampled values X(m) producing the
% time-domain sequence x(n).
%
% If the sequence, X(m), length is less than N, then the sequence is padded
% with zeros.  If the sequence, X(m), is greater than N, then the sequence
% is truncated to N elements.
%
% This algorithm uses the (forward) radix-2 FFT to calculate the inverse.
%
% @param   X      The discrete sequence of frequency-domain sampled values.
% @param   N      The length of the FFT.
%
% @return  The time-domain sequence x(n).
%
% @throw   Error if the input sequence is not a vector.
% @throw   Error if N is not a power of two.
%
% [1] Lyons, Richard G., Understanding Digital Signal Processing 2nd 
%     Edition, Prentice Hall 2004
% [2] http://en.wikipedia.org/wiki/Cooley?Tukey_FFT_algorithm
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
function [x] = ir2fft(X,N)

    if ~isvector(X)
        error('ir2fft: input sequence must be a vector');
    end

    if ~exist('N','var')
        N = nextpow2(length(X));
    end
 
    % ensure that N is a power of 2
    N  = abs(ceil(N(1)));
    Np = 2^nextpow2(N);

    if N~=Np
        fprintf(['ir2fft: FFT length %d is not a power of 2; '...
                 'converting to %d\n'],N,Np);
        N = Np;
    end

    % do we need to pad x or truncate it?
    ctype = class(X);
    diff  = N-length(X);
    if diff > 0
        if isrow(X)
            X = [X cast(zeros(1,diff),ctype)];
        else
            X = [X;cast(zeros(diff,1),ctype)];
        end
    elseif diff < 0
        diff = abs(diff);
        X = X(1:end-diff);
    end
   
    x = conj(r2fft(conj(X),N))./N;
     
end
