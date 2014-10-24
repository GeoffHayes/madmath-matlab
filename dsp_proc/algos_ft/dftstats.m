%**************************************************************************
% Returns the magnitude, power, and phase for each value in the
% DFT-generated frequency domain sequence.
%
% @param   X   The frequency-domain sequence generated from the DFT.
% @param   fs  The sampling frequency (samples per second) that the data
%              used to generate X was sampled at.
% @param   x   The time-domain sequence that generated X via the DFT.
%
% @return  Statistics for each value in the frequency-domain sequence: the
%          magnitude, power, phase angle (degrees), frequencies (Hz).
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
function [STATS] = dftstats(X,fs,x)

    if ~isvector(X)
        error('dftstats: input sequence must be a vector');
    end

    realX = real(X);
    imagX = imag(X);
    
    realX(abs(realX)<0.000000001)=0;
    imagX(abs(imagX)<0.000000001)=0;

    N = length(X);
    
    tdIsReal  = 0;
    ignoreClr = [255 224 224]./255;
    if exist('x','var')
        tdIsReal = isreal(x);
    end
    
    if tdIsReal
        % if the time-domain sequence data is real, then the magnitude
        % factor of the DFT data is:
        magFactor = N/2;
    else
        % else the time-domain sequence data is complex and so the
        % magnitude factor of the DFT data is:
        magFactor = N;
    end
        
    
    STATS = zeros(N,4);
    
    % compute the magnitude for each value
    STATS(:,1) = abs(X);
    
    % compute the power for each value
    STATS(:,2) = STATS(:,1).^2;
    
    % compute the phase angle for each value
    STATS(:,3) = atan2d(imagX,realX);
    
    % compute the frequencies for each value
    STATS(:,4) = 0:fs/N:(N-1)*fs/N;
    
    figure;
    
    % plot the magnitude
    subplot(2,2,1);
    hold on;
    freq = STATS(:,4);
    
    % if the time-domain sequence is real then we can ignore half of the
    % output of the DFT...so colour this region starting at N/2+2 as:
    % 0<=m<=N/2 --> 1<=m<=N/2+1 (in one-based arrays) as the valid
    % non-duplicated data
    if tdIsReal
        fill([freq(ceil(N/2)+2); freq(end)+fs/N; freq(end)+fs/N    ; freq(ceil(N/2)+2)], ...
             [0                ; 0             ; max(STATS(:,1))   ; max(STATS(:,1))  ], ...
             ignoreClr);
    end    
    
    plot(freq,STATS(:,1),'b.','MarkerSize',8);
    title('Magnitude of X(m)');
    xlabel('Frequency (Hz)');
    for i=1:N
        line([freq(i) freq(i)],[0 STATS(i,1)],'Color','b'); 

        if floor((STATS(i,1)+0.5)/magFactor) > 0
            % plot the true amplitude of the time-domain signal from their
            % DFT spectral result (only if non-zero)
            plot(freq(i),STATS(i,1)/magFactor,'r.','MarkerSize',8);
        end
       
    end
    line([freq(1) freq(end)+fs/N],[0 0],'Color','k');
    axis([-1 freq(end)+fs/N 0 max(STATS(:,1))]);

    % plot the real part of X(m)
    subplot(2,2,2);
    hold on;
    
    % if the time-domain sequence is real then we can ignore half of the
    % output of the DFT...so colour this region starting at N/2+2 as:
    % 0<=m<=N/2 --> 1<=m<=N/2+1 (in one-based arrays) as the valid
    % non-duplicated data
    if tdIsReal
        fill([freq(ceil(N/2)+2); freq(end)+fs/N; freq(end)+fs/N    ; freq(ceil(N/2)+2)], ...
             [min(realX)       ; min(realX)    ; max(realX)+eps    ; max(realX)+eps   ], ...
             ignoreClr);
    end    
    
    plot(freq,realX,'b.','MarkerSize',8);
    title('Real part of X(m)');
    xlabel('Frequency (Hz)');
    for i=1:N
       line([freq(i) freq(i)],[0 realX(i)],'Color','b'); 
    end
    line([freq(1) freq(end)+fs/N],[0 0],'Color','k');
    axis([-1 freq(end)+fs/N min(realX) max(realX)+eps]);
    
    % plot the phase
    subplot(2,2,3);
    hold on;
    freq = 0:fs/N:(N-1)*fs/N;
    
    % if the time-domain sequence is real then we can ignore half of the
    % output of the DFT...so colour this region starting at N/2+2 as:
    % 0<=m<=N/2 --> 1<=m<=N/2+1 (in one-based arrays) as the valid
    % non-duplicated data
    if tdIsReal
        fill([freq(ceil(N/2)+2); freq(end)+fs/N; freq(end)+fs/N    ; freq(ceil(N/2)+2)], ...
             [-180             ; -180          ; 180               ; 180              ], ...
             ignoreClr);
    end     
    
    plot(freq,STATS(:,3),'b.','MarkerSize',8);
    title('Phase of X(m) (degs)');
    xlabel('Frequency (Hz)');
    for i=1:N
       line([freq(i) freq(i)],[0 STATS(i,3)],'Color','b'); 
    end
    line([freq(1) freq(end)+fs/N],[0 0],'Color','k');
    axis([-1 freq(end)+fs/N -180 180]);  
    
    % plot the imaginary part of X(m)
    subplot(2,2,4);
    hold on;
    
    % if the time-domain sequence is real then we can ignore half of the
    % output of the DFT...so colour this region starting at N/2+2 as:
    % 0<=m<=N/2 --> 1<=m<=N/2+1 (in one-based arrays) as the valid
    % non-duplicated data
    if tdIsReal
        fill([freq(ceil(N/2)+2); freq(end)+fs/N; freq(end)+fs/N    ; freq(ceil(N/2)+2)], ...
             [min(imagX)       ; min(imagX)    ; max(imagX)+eps    ; max(imagX)+eps   ], ...
             ignoreClr);
    end       
    
    plot(freq,imagX,'b.','MarkerSize',8);
    title('Imaginary part of X(m)');
    xlabel('Frequency (Hz)');
    for i=1:N
       line([freq(i) freq(i)],[0 imagX(i)],'Color','b'); 
    end
    line([freq(1) freq(end)+fs/N],[0 0],'Color','k');
    axis([-1 freq(end)+fs/N min(imagX) max(imagX)+eps]);
   
end

% from [1], section 3.11: DFT frequency domain sampling - increasing the
% block size by padding (with zeros) the sampled data lets us better model
% or approximate the continuous fourier transform of the function (x1)
% which is sampled for the time period T only and zeros elsewhere:

% x1 = @(t)(sin(2*pi*3*t));
% close all;N=16;fs=16;[A,D]=wava2d(fs,fs,x1);X=dft(D,N);STATS=dftstats(X,N);
% close all;N=32;fs=16;[A,D]=wava2d(fs,fs,x1);X=dft(D,N);STATS=dftstats(X,N);
% close all;N=64;fs=16;[A,D]=wava2d(fs,fs,x1);X=dft(D,N);STATS=dftstats(X,N);
% close all;N=128;fs=16;[A,D]=wava2d(fs,fs,x1);X=dft(D,N);STATS=dftstats(X,N);

% while the DFT frequency-domain sampling characteristic becomes more
% obvious, the bin index for the centre of the main lobe is different for
% each of the DFT outputs.  the centre frequency for the mth bin is equal
% to m*fs/N

