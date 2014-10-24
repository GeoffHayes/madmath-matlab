%**************************************************************************
% Performs the radix-2 Fast Fourier Transform (FFT) on the discrete 
% sequence of time-domain sampled values x[n] of a continuous variable x(y), 
% producing the frequency-domain sequence X(m).
%
% If the sequence, x(n), length is less than N, then the sequence is padded
% with zeros.  If the sequence, x(n), is greater than N, then the sequence
% is truncated to N elements.
%
% This algorithm is known as decimation in time since the input sequence is
% continually divided into its even and odd parts.
%
% @param   x      The discrete sequence of time-domain sampled values.
% @param   N      The length of the FFT.
%
% @return  The frequency-domain sequence X(m).
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
function [X] = r2fft(x,N)

    % frequency domain data to be populated via the recursive function
    global fdData;

    if ~isvector(x)
        error('r2fft: input sequence must be a vector');
    end

    if ~exist('N','var')
        N = length(x);
    end
 
    % ensure that N is a power of 2
    N  = abs(ceil(N(1)));
    Np = 2^nextpow2(N);

    if N~=Np
        fprintf(['r2fft: FFT length %d is not a power of 2; '...
                 'converting to %d\n'],N,Np);
        N = Np;
    end
    
    piFctrs = 2*pi*(0:1:N-1)';

    % do we need to pad x or truncate it?
    ctype = class(x);
    diff  = N-length(x);
    if diff > 0
        if isrow(x)
            x = [x cast(zeros(1,diff),ctype)];
        else
            x = [x;cast(zeros(diff,1),ctype)];
        end
    elseif diff < 0
        diff = abs(diff);
        x = x(1:end-diff);
    end
    
    % size the frequency domain data
    fdData = zeros(N,1);
   
    % perform the FFT
    r2fft_recurse(x,piFctrs,N,1,1,1);
    
    % set and return the data
    X = fdData;
     
end

% recursive helper function for the radix-2 FFT
function r2fft_recurse(x,piFctrs,N,start,step,u)

    global fdData;
    
    if N==2
        fdData(u)   = x(start);
        fdData(u+1) = x(start+step);
    else
        r2fft_recurse(x,piFctrs,N/2,start,2*step,u);
        r2fft_recurse(x,piFctrs,N/2,start+step,2*step,N/2+u);
    end

    for k=0:N/2-1
        t        = fdData(u+k);
        a        = piFctrs(1+k)/N;
        b        = (cos(a) - (1i)*sin(a))*fdData(u+k+N/2);
        fdData(u+k)     = t + b;
        fdData(u+k+N/2) = t - b;
    end 
end
