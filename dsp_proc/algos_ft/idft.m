%**************************************************************************
% Performs an Inverse Discrete Fourier Transform (IDFT) on the frequency
% domain sequence X(m) output from the DFT, producing the time-domain
% sequence x(n) where:
%
%           x(n) = (1/N)SUM{m=0...N-1}X(m)exp(i2pimn/N)
%
% If the sequence, X(m), length is less than N, then the sequence is padded
% with zeros.  If the sequence, X(m), is greater than N, then the sequence
% is truncated to N elements.
%
% where N is the length of the sequence x.
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
function [x] = idft(X,N)

    if ~isvector(X)
        error('idft: input sequence must be a vector');
    end
    
    if ~exist('N','var')
        N = length(X);
    else
        N = abs(ceil(N(1)));
        
        % do we need to pad X or truncate it?
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
    end
    
    % from [1], with N input frequency-domain sample values, the IDFT
    % reverses the DFT process which produced the spectral content of the 
    % time-domain input at N equally spaced frequency points.
    %
    
    % using Euler's relationship, exp(-ja) = cos(a)-jsin(a) the DFT
    % equation can be rewritten as:
    %
    %           x(n) = SUM{m=0...N-1}X(m)[cos(2pinm/N)+isin(2pinm/N)]
    %
    x = zeros(N,1);
    
    for n=0:N-1
        sum = 0;
        for m=0:N-1
            sum = sum + X(m+1)*(cos(2*pi*m*n/N) + 1i*sin(2*pi*m*n/N));
        end
        
        x(n+1) = sum/N;
    end
end

% as with the DFT, linearity (and other properties) holds true [1]

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

% x=idft(X);

