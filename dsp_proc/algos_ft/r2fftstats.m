%**************************************************************************
% Returns the magnitude, power, and phase for each value in the
% FFT-generated frequency domain sequence.
%
% @param   X         The frequency-domain sequence generated from the FFT.
% @param   fs        The sampling frequency (samples per second) that the data
%                    used to generate X was sampled at.
% @param   x         The time-domain sequence that generated X via the FFT.
% @param   hs        A vector of axes handles to draw the data on (optional).
% @param   maxFrqHz  The maximum frequency (Hz) shown in each display.
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
function [STATS] = r2fftstats(X,fs,x,hs,maxFrqHz)

    if ~isvector(X)
        error('r2fftstats: input sequence must be a vector');
    end
    
    if ~exist('hs','var')
        hs = [];
    end
    
    if ~exist('maxFrqHz','var')
        maxFrqHz = fs;
    else
        maxFrqHz = min(fs,maxFrqHz);
    end

    realX = real(X);
    imagX = imag(X);
    
    realX(abs(realX)<0.000000001)=0;
    imagX(abs(imagX)<0.000000001)=0;

    N = length(X);
    
    tdIsReal  = 0;
    ignoreClr = [128 128 128]./255;
    if exist('x','var')
        tdIsReal = isreal(x);
    end
    
    if tdIsReal
        % if the time-domain sequence data is real, then the magnitude
        % factor of the FFT data is:
        magFactor = N/2;
    else
        % else the time-domain sequence data is complex and so the
        % magnitude factor of the FFT data is:
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

    if isempty(hs)
        figure;
        h1 = subplot(2,2,1);
    else
        h1 = hs(1);
    end
    
    freq = STATS(:,4);
    
    idx  = find(freq>maxFrqHz,1)-1;
    if isempty(idx)
        idx = length(freq);
    end
    freq = freq(1:idx);
    
    if h1>0
    
        cla(h1);

        % plot the magnitude

        % if the time-domain sequence is real then we can ignore half of the
        % output of the FFT...so colour this region starting at N/2+2 as:
        % 0<=m<=N/2 --> 1<=m<=N/2+1 (in one-based arrays) as the valid
        % non-duplicated data
        if tdIsReal && length(freq)>(ceil(N/2)+2)
            fill([freq(ceil(N/2)+2); freq(end)+fs/N; freq(end)+fs/N    ; freq(ceil(N/2)+2)], ...
                 [0                ; 0             ; max(STATS(:,1))   ; max(STATS(:,1))  ], ...
                 ignoreClr, 'Parent',h1);
        end    

        hold(h1,'on');

        stem(h1,freq,STATS(1:idx,1),'b','MarkerSize',3,'MarkerFaceColor','b');
    %     idcs = find(floor(STATS(:,1)+0.5/magFactor)>0);
    %     if ~isempty(idcs)
    %         try
    %         hold(h1,'on');
    %         scatter(h1,freq(idcs),STATS(idcs,1)/magFactor,'r');
    %         catch exc
    %             t=3;
    %         end
    %     end
        title(h1,'Magnitude of X(m)');
        xlabel(h1,'Frequency (Hz)');
        %line([freq(1) freq(end)+fs/N],[0 0],'Color','k','Parent',h1);
        axis(h1,[-1 freq(end)+fs/N 0 max(STATS(:,1))+eps]);
        
    end

    % plot the real part of X(m)
    if isempty(hs)
        h2 = subplot(2,2,2);
    elseif length(hs)>1
        h2 = hs(2);
    else
        h2 = [];
    end

    % if the time-domain sequence is real then we can ignore half of the
    % output of the FFT...so colour this region starting at N/2+2 as:
    % 0<=m<=N/2 --> 1<=m<=N/2+1 (in one-based arrays) as the valid
    % non-duplicated data
    if ~isempty(h2) && h2>0
        
        cla(h2);
        
        maxAbsRealX = max(abs(min(realX)),abs(max(realX)));
        
        if tdIsReal && length(freq)>(ceil(N/2)+2)
            fill([freq(ceil(N/2)+2); freq(end)+fs/N; freq(end)+fs/N    ; freq(ceil(N/2)+2)], ...
                 [-maxAbsRealX       ; -maxAbsRealX   ; maxAbsRealX    ; maxAbsRealX      ], ...
                 ignoreClr, 'Parent',h2);
        end   
        
        hold(h2,'on');
        
        stem(h2,freq,realX(1:idx),'b','MarkerSize',1);
        title(h2,'Real part of X(m)');
        xlabel(h2,'Frequency (Hz)');
        %line([freq(1) freq(end)+fs/N],[0 0],'Color','k','Parent',h2);
        axis(h2,[-1 freq(end)+fs/N -maxAbsRealX maxAbsRealX+eps]);
    end
    
    % plot the phase
    if isempty(hs)
        h3 = subplot(2,2,3);
    elseif length(hs)>2
        h3 = hs(3);
    else
        h3 = [];
    end
    
    freq = 0:fs/N:(N-1)*fs/N;
    idx  = find(freq>maxFrqHz,1)-1;
    if isempty(idx)
        idx = length(freq);
    end
    freq = freq(1:idx);

    if ~isempty(h3) && h3>0
        
        cla(h3);
        
        % if the time-domain sequence is real then we can ignore half of the
        % output of the FFT...so colour this region starting at N/2+2 as:
        % 0<=m<=N/2 --> 1<=m<=N/2+1 (in one-based arrays) as the valid
        % non-duplicated data
        if tdIsReal && length(freq)>(ceil(N/2)+2)
            fill([freq(ceil(N/2)+2); freq(end)+fs/N; freq(end)+fs/N    ; freq(ceil(N/2)+2)], ...
                 [-180             ; -180          ; 180               ; 180              ], ...
                 ignoreClr, 'Parent', h3);
        end  
        
        hold(h3,'on');
        
        stem(h3,freq,STATS(1:idx,3),'b','MarkerSize',1);
        title(h3,'Phase of X(m) (degs)');
        xlabel(h3,'Frequency (Hz)');
        %line([freq(1) freq(end)+fs/N],[0 0],'Color','k','Parent',h3);
        axis(h3,[-1 freq(end)+fs/N -180 180]);  
    end
    
    % plot the imaginary part of X(m)
    if isempty(hs)
        h4 = subplot(2,2,4);
    elseif length(hs)>3
        h4 = hs(4);
    else
        h4 = [];
    end
    
    if ~isempty(h4) && h4>0
        
        cla(h4);
        
        maxAbsImagX = max(abs(min(imagX)),abs(max(imagX)));
        
        % if the time-domain sequence is real then we can ignore half of the
        % output of the FFT...so colour this region starting at N/2+2 as:
        % 0<=m<=N/2 --> 1<=m<=N/2+1 (in one-based arrays) as the valid
        % non-duplicated data
        if tdIsReal && length(freq)>(ceil(N/2)+2)
            fill([freq(ceil(N/2)+2); freq(end)+fs/N; freq(end)+fs/N    ; freq(ceil(N/2)+2)], ...
                 [-maxAbsImagX     ; -maxAbsImagX   ; maxAbsImagX   ; maxAbsImagX         ], ...
                 ignoreClr, 'Parent', h4);
        end  
        
        hold(h4,'on');
        
        stem(h4,freq,imagX(1:idx),'b','MarkerSize',1);
        title(h4,'Imaginary part of X(m)');
        xlabel(h4,'Frequency (Hz)');
        %line([freq(1) freq(end)+fs/N],[0 0],'Color','k','Parent',h4);
        
        axis(h4,[-1 freq(end)+fs/N -maxAbsImagX maxAbsImagX+eps]);
    end
   
end

% from [1], section 3.11: FFT frequency domain sampling - increasing the
% block size by padding (with zeros) the sampled data lets us better model
% or approximate the continuous fourier transform of the function (x1)
% which is sampled for the time period T only and zeros elsewhere:

% x1 = @(t)(sin(2*pi*3*t));
% close all;N=16;fs=16;[A,D]=wava2d(fs,fs,x1);X=r2fft(D,N);STATS=r2fftstats(X,fs);
% close all;N=32;fs=16;[A,D]=wava2d(fs,fs,x1);X=r2fft(D,N);STATS=r2fftstats(X,fs);
% close all;N=64;fs=16;[A,D]=wava2d(fs,fs,x1);X=r2fft(D,N);STATS=r2fftstats(X,fs);
% close all;N=128;fs=16;[A,D]=wava2d(fs,fs,x1);X=r2fft(D,N);STATS=r2fftstats(X,fs);

% while the FFT frequency-domain sampling characteristic becomes more
% obvious, the bin index for the centre of the main lobe is different for
% each of the FFT outputs.  the centre frequency for the mth bin is equal
% to m*fs/N

