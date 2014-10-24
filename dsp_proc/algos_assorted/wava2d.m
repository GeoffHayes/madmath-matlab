%**************************************************************************
% Plots a discrete time signal of a wave given an input wave function and
% sampling period.  This simulates the passing of the analog wave through
% an analog-to-digital (A/D) converter.  This results in the analog signal
% being sampled.
%
% @param   N         The size of the block to create (can be the DFT/FFT 
%                    size).
% @param   fs        The sampling rate in samples per second.
% @param   wavFunc   The wave function.
% @param   wdwFunc   The optional window function.
%
% @return  A         The simulated analog data.
% @return  D         The digital/discrete data sampled from the analog
%                    "input".
%
% @throw   Error if the block size (blkSize) is not a scalar.
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
function [A,D] = wava2d(N,fs,wavFunc,wdwFunc)

    if ~isscalar(N)
        error('wava2d: block size must be a scalar');
    end
    
    if ~exist('wdwFunc','var')
        wdwFunc = [];
    end
    
    fs = abs(fs);
    N  = abs(ceil(N));
    
    % use a gap of 1/1000th of a second
    analogGap = 0.0001;

    % if the sampling rate is too great, fix to the max allowed (which is the
    % gap)
    numSamplesPerSecond = 1/analogGap;
    if fs > numSamplesPerSecond
        fs = numSamplesPerSecond;
    end    
    ts = 1/fs;

    Td = (0:ts:N*ts-ts)';
    Tc = (0:analogGap:Td(end))';

    % given the time, compute the input analog time domain data evaluated
    % for the input function
    anaTimeDomainData = wavFunc(Tc);
    
    % apply the A/D converter against the analog time domain data
    digTimeDomainData = wavFunc(Td);
    
    A = anaTimeDomainData;
    D = digTimeDomainData;

    figure;
    
    numYPlots = 3;
    numXPlots = 1;
    
    if ~isempty(wdwFunc)
        numYPlots = numYPlots + 1;
    end
    
    dataIsComplex = 0;
    if ~isreal(D)
        numXPlots     = numXPlots + 1;
        dataIsComplex = 1;
    end
    
    maxNumPlots = numXPlots*numYPlots;
    
    atPlot = 1;

    % plot the continuous wave
    subplot(numYPlots,numXPlots,atPlot);
    hold on;
    line([Tc(1) Tc(end)],[0 0],'Color','k');
    xlabel('Time (sec)');    
    atPlot = atPlot + 1;
    if dataIsComplex
        rA = real(A);
        title('Analog Data (real)');
        plot(Tc,rA);
        axis([0 Tc(end)+ts/10 min(rA)-0.5 max(rA)+0.5]);
    else
        title('Analog Data');
        plot(Tc,A);
        axis([0 Tc(end)+ts/10 min(A)-0.5 max(A)+0.5]);
    end
    
    if dataIsComplex
        subplot(numYPlots,numXPlots,atPlot);
        hold on;
        atPlot = atPlot + 1;
        iA = imag(A);
        title('Analog Data (imag)');
        plot(Tc,iA);
        line([Tc(1) Tc(end)],[0 0],'Color','k');
        xlabel('Time (sec)');
        axis([0 Tc(end)+ts/10 min(iA)-0.5 max(iA)+0.5]);  
    end

    % note that if the frequency, f0, is equal to f0+k*fs, then the discrete
    % data is indistinguishable for any positive or negative integer k.  Thus
    % the sampled sinusoid can represent an infinite number of (continuous)
    % sinusoids.  this is because:

    %   sin(2*pi*(f0+k*fs)*n*ts) = sin(2*pi*f0*n*ts + 2*pi*k*fs*n*ts)
    % = sin(2*pi*f0*n*ts)*cos(2*pi*k*fs*n*ts) + sin(2*pi*k*fs*n*ts)*cos(2*pi*f0*n*ts)
    % = sin(2*pi*f0*n*ts)*1 + 0*cos(2*pi*f0*n*ts)
    % = sin(2*pi*f0*n*ts)
    % = x(n)

    % in the above, we could compute the discreteData as:
    % sinee(2*pi*discreteTime*frequency) but since we already have the data..

    % plot the discrete-time signal given the sampling rate
    subplot(numYPlots,numXPlots,atPlot);
    hold on;
    line([Td(1) Td(end)],[0 0],'Color','k');
    xlabel('Time (sec)');
    atPlot = atPlot + 1;
    if dataIsComplex
        rD = real(D);
        title(['Digital Data (real) at ' num2str(fs) ...
            ' samples per second']);
        plot(Td,rD,'r.', 'MarkerSize',8);
        axis([0 Tc(end)+ts/10 min(rD)-0.5 max(rD)+0.5]);
        
        % draw a lollipop plot
        for i=1:length(Td)
            line([Td(i) Td(i)],[0 rD(i)],'Color','r'); 
        end        
    else
        title(['Digital Data at ' num2str(fs) ...
            ' samples per second']);
        plot(Td,D,'r.', 'MarkerSize',8);
        axis([0 Tc(end)+ts/10 min(D)-0.5 max(D)+0.5]);
        
        % draw a lollipop plot
        for i=1:length(Td)
            line([Td(i) Td(i)],[0 D(i)],'Color','r'); 
        end
    end    
    
    if dataIsComplex
        subplot(numYPlots,numXPlots,atPlot);
        hold on;
        atPlot = atPlot + 1;
        iD = imag(D);
        title(['Digital Data (imag) at ' num2str(fs) ...
            ' samples per second']);
        plot(Td,iD,'r.', 'MarkerSize',8);
        line([Td(1) Td(end)],[0 0],'Color','k');
        xlabel('Time (sec)');
        axis([0 Tc(end)+ts/10 min(iD)-0.5 max(iD)+0.5]);  
        
        % draw a lollipop plot
        for i=1:length(Td)
            line([Td(i) Td(i)],[0 iD(i)],'Color','r'); 
        end         
    end   
    
    % draw the windowed plots if needed
    if numYPlots==4
        subplot(numYPlots,numXPlots,atPlot);
        hold on;
        line([Td(1) Td(end)],[0 0],'Color','k');
        xlabel('Time (sec)');
        atPlot = atPlot + 1;

        % compute the window weights and apply them against the discrete
        % data
        w = wdwFunc(N);
        D = w.*D;
        
        if dataIsComplex
            rD = real(D);
            title(['Windowed Digital Data (real) at ' num2str(fs) ...
                ' samples per second']);
            plot(Td,rD,'r.', 'MarkerSize',8);
            axis([0 Tc(end)+ts/10 min(rD)-0.5 max(rD)+0.5]);

            % draw a lollipop plot
            for i=1:length(Td)
                line([Td(i) Td(i)],[0 rD(i)],'Color','r'); 
            end        
        else
            title(['Windowed Digital Data at ' num2str(fs) ...
                ' samples per second']);
            plot(Td,D,'r.', 'MarkerSize',8);
            axis([0 Tc(end)+ts/10 min(D)-0.5 max(D)+0.5]);

            % draw a lollipop plot
            for i=1:length(Td)
                line([Td(i) Td(i)],[0 D(i)],'Color','r'); 
            end
        end   
        
        % draw the window
        plot(Td,w,'b-');        

        if dataIsComplex
            subplot(numYPlots,numXPlots,atPlot);
            hold on;
            atPlot = atPlot + 1;
            iD = imag(D);
            title(['Windowed Digital Data (imag) at ' num2str(fs) ...
                ' samples per second']);
            plot(Td,iD,'r.', 'MarkerSize',8);
            line([Td(1) Td(end)],[0 0],'Color','k');
            xlabel('Time (sec)');
            axis([0 Tc(end)+ts/10 min(iD)-0.5 max(iD)+0.5]);  

            % draw a lollipop plot
            for i=1:length(Td)
                line([Td(i) Td(i)],[0 iD(i)],'Color','r'); 
            end    
            
            % draw the window
            plot(Td,w,'b-');
        end        
    end    

    % now do an FFT on the discrete data
    L    = length(D);
    NFFT = 2^nextpow2(L);
    Y    = fft(D,NFFT)/L;
    f    = fs/2*linspace(0,1,NFFT/2+1);
    
    subplot(numYPlots,numXPlots,atPlot:maxNumPlots);
    plot(f,2*abs(Y(1:NFFT/2+1))) 
    title('Single-Sided Amplitude Spectrum of y(t)')
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')

end

% usage - define a sine wave for some frequency and add (or not) some
% gaussian noise

% x1 = @(t)(sin(2*pi*1*t)+2*rand(size(t)));
% [A,D] = wava2d(512,100,x1);

% x1 = @(t)(0.7*sin(2*pi*50*t) + sin(2*pi*120*t) + 2*rand(size(t)));
% fs = 2*120+5;
% [A,D] = wava2d(512,fs,x1);

% Note that the 2kHz term is shifted in phase by 35 degrees relative to the
% 1kHz sinewave:
% x1 = @(t)(sin(2*pi*1000*t) + 0.5*sin(2*pi*2000*t+3*pi/4));

% complex usage:
% x2 = @(t)exp((1i)*2*pi*2*t);
% [A,D] = wava2d(64,64,x2);
% [A,D] = wava2d(64,64,x2,@wdw_triangle);