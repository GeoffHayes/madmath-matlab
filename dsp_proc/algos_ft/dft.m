%**************************************************************************
% Performs a Discrete Fourier Transform (DFT) on the discrete sequence of
% time-domain sampled vales x[n] of a continuous variable x(y), producing
% the frequency-domain sequence X(m) where:
%
%                  X(m) = SUM{n=0...N-1}x(n)exp(-i2pinm/N)
%
% If the sequence, x(n), length is less than N, then the sequence is padded
% with zeros.  If the sequence, x(n), is greater than N, then the sequence
% is truncated to N elements.
%
% @param   x      The discrete sequence of time-domain sampled values.
% @param   N      The length of the DFT.
%
% @return  The frequency-domain sequence.
%
% @throw   Error if the input sequence is not a vector.
%
% [1] Lyons, Richard G., Understanding Digital Signal Processing 2nd 
%     Edition, Prentice Hall 2004
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
function [X] = dft(x,N,shift)

    if ~isvector(x)
        error('dft: input sequence must be a vector');
    end
    
    if ~exist('N','var')
        N = length(x);
    else
        N = abs(ceil(N(1)));
        
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
    end
    
    if ~exist('shift','var')
        shift = 0;
    elseif ~isscalar(shift)
        error('dft: shift input is not a scalar');
    else
        shift = floor(shift);
    end
    
    % from [1], with N input time-domain sample values, the DFT determines
    % the spectral content of the input at N equally spaced frequency 
    % points.  The value N is an important parameter because it determines
    % how many input samples are needed, then resolution of the
    % frequency-domain results, and the amount of processing time necessary
    % to calculate an N-point DFT
    %
    % from [1], the j (sqrt(-1)) is an abstraction that helps us compare
    % the phase relationship between various sinusoidal components of a
    % signal
    
    % using Euler's relationship, exp(-ja) = cos(a)-jsin(a) the DFT
    % equation can be rewritten as:
    %
    %           X(m) = SUM{n=0...N-1}x(n)[cos(2pinm/N)-isin(2pinm/N)]
    %
    X = zeros(N,1);
    
    for m=0:N-1
        sum = 0;
        for n=0:N-1
            sum = sum + x(n+1)*(cos(2*pi*n*m/N) - 1i*sin(2*pi*n*m/N));
        end
        % apply the shift
        X(m+1) = sum*exp(1i*2*pi*shift*m/N);
        
        % from [1], each X(m) DFT output term is the sum of the point for
        % point product between an input sequence of signal values and a
        % complex sinusoid of th form cos(a)-jsin(a).  The exact
        % frequencies of the different sinusoids depend on both the
        % sampling rate fs at which the original signal was sampled and the
        % number of samples N.  For example, if we are sampling a
        % continuous signal at a rate of 500 samples/sec and then perform a
        % N=16-point DFT on the sampled data, the fundamental frequcy of
        % the sinusoids is fs/N=500/16 or 31.25 Hz.  The other X(m)
        % analysis frquencies are integral multiples of the fundamental
        % frequency:
        % X(0) = first frequency term with analysis frequency 0*31.25 = 0Hz
        % X(1) = second .....                                 1*31.25 = 31.25Hz
        % X(2) = third  .....                                 2*31.25 = 62.50Hz
        % etc.
        %
        % The N separate DFT analysis frequencies are m*fs/N for m=0...N-1
        %
        % The X(0) DFT term tells us the magnitude of any 0-Hz (DC)
        % component in the input signal.
        % The X(k) DFT term tells us the magnitude of any k*31.25Hz
        % component in the input signal.
    end
    
    % note that if x is a vector of real inputs only, then ([1]) the
    % complex outputs for m=1...N/2-1 are redundant with frequency output
    % values for m>(N/2).  The mth DFT output will ahve the same magnitude
    % as the (N-M)th DFT ouptut.  The phase angle of the DFT's mth output
    % is the negative of the phase angle of the (N-m)th DFT output:
    
    % X(m)   = SUM{n=0...N-1}x(n)exp(-i2pinm/N)
    % X(N-m) = SUM{n=0...N-1}x(n)exp(-i2pin(N-m)/N)
    %        = SUM(n=0...N-1}x(n)exp(-i2pin(N)/N)exp(-i2pin(-m)/N)
    %        = SUM(n=0...N-1}x(n)exp(-i2pin)exp(i2pinm/N)
    % note that exp(-i2pin) = cos(2pin)-isin(2pin) = 1 for all n
    % X(N-m) = SUM(n=0...N-1}x(n)exp(i2pinm/N)
    % with the only difference from X(m) being the lack of the negative on
    % the exp product.  This is just the complex conjugate of X(m) which
    % explains the behaviour of the negative phases and negative imaginary
    % parts in the third and fourth suplots.
    
    % again from [1], the DFT has the linearity property in that:
    %           DFT(x1(n) + x2(n)) = DFT(x1(n)) + DFT(x2(n))
    % This can be proven given the exponential DFT equation.
    
    % Note that if the input to the DFT is real, the output magnitude of
    % the DFT, for the sinewave at a certain frequency with an amplitude of
    % A, is A*N/2.  If the input is complex, then the output magnitude of
    % the DFT is A*N.
    
    % the frequency resolution of the DFT is fs/N.  This will indicate the
    % frequency associated with the large DFT magnitude term
    
end

% usage: the following sinusoidal signal has three frequencies at 1000Hz,
% 2000Hz and 125Hz

% x1 = @(t)(sin(2*pi*1000*t) + 0.5*sin(2*pi*2000*t+3*pi/4) + 0.75*sin(2*pi*125*t));
% x1 = @(t)(sin(2*pi*1000*t) + 0.5*sin(2*pi*2000*t+3*pi/4) + 0.75*sin(2*pi*125*t) + 0.0*rand(size(t)));

% in order to find the 1000Hz signal the FFT block size (N) must be such 
% that there exists a k (k=0,1,2,...) with fs/N*k = 1000.  If fs=8000, then
% N=8,k=1 or N=16,k=2 or N=32,k=4 etc.

% to find all three frequencies, then N must be chosen such that there
% exists a u,v,w with fs/N*u=1000, fs/N*v=2000, and fs/N*w=125 i.e.
% N=64,128, etc.

% close all;[A,D]=wava2d(8,8000,x1);X=dft(D);STATS=dftstats(X,8000);

% note that if the sample rate equals the block size, then whatever
% frequency the sine wave is operating at is what will be sampled.  So a
% 3Hz sine wave with N=fs=64 will have all three cycles sampled.  If
% N=2*fs, then six cycles of the sine wave will be sampled.  Etc.

% DFT leakage example because the input sequence x1 does not have an
% integral number of cycles over our 64-sample interval, and so energy has
% leaked into all the other DFT output bins:

% x1 = @(t)(sin(2*pi*3.4*t));

% close all;N=64,fs=64;[A,D]=wava2d(N,fs,x1);X=dft(D,N);STATS=dftstats(X,fs);

