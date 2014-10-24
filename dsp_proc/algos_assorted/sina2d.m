%**************************************************************************
% Plots a discrete time signal of a sine wave given a frequency and period.
% This simulates the passing of the analog sine wave through an
% analog-to-digital (A/D) converter.  This results in the analog signal
% being sampled.
%
% @param   freqHz    The frequency of the sine wave in hertz (number of
%                    cycles of the sine wave per second).
% @param   plotSec   The period in seconds to plot the data.
% @param   fs        The sampling rate in samples per second.
%
% @return  anaData   The simulated analog data.
% @return  digData   The digital/discrete data sampled from the analog
%                    "input".
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
function [anaData,digData] = sina2d(freqHz,plotSec,fs)

    freqHz  = abs(freqHz);
    plotSec = abs(plotSec);
    fs      = abs(fs);

    % use a gap of 1/1,000,000th of a second
    analogGap = 0.000001;

    anaTime = 0:analogGap:plotSec;
    anaData = sin(2*pi*anaTime*freqHz);

    figure;

    % plot the continuous sine wave
    subplot(2,1,1);
    plot(anaTime,anaData);
    line([anaTime(1) anaTime(end)],[0 0],'Color','k');
    title(['Analog Data (' num2str(freqHz) ' Hz)']);
    xlabel('Time (sec)');
    anaXTicks = get(gca,'XTickLabel');

    % if the sampling rate is too great, fix to the max allowed (which is the
    % gap)
    numSamplesPerSecond = 1/analogGap;
    if fs > numSamplesPerSecond
        fs = numSamplesPerSecond;
    end

    ts = 1/fs;

    % apply the A/D converter against the sine wave (analog data)
    digTime = 0:ts:plotSec;                  % n
    digData = sin(2*pi*freqHz*digTime);      % x(n)

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
    subplot(2,1,2);
    plot(digTime,digData,'r.', 'MarkerSize',8);
    line([digTime(1) digTime(end)],[0 0],'Color','k');
    title(['Digital/Discrete Data (' num2str(freqHz) ' Hz) at ' num2str(fs) ...
        ' samples per second']);
    xlabel('Time (sec)');
    
    % draw a lollipop plot
    for i=1:length(digTime)
       line([digTime(i) digTime(i)],[0 digData(i)],'Color','r'); 
    end
    
    set(gca,'XTickLabel',anaXTicks);

    % note that the sampling rate, fs, is in samples per second
    % this means that the sample period, ts, is 1/fs in seconds per sample

    % if the number of samples per period (i.e. those samples that belong to a
    % single period before they start repeating is x then we can estimate the
    % sinewave period via:
    %
    % sinewave period = x samples/period * ts
    %
    % The units for the sinewave period is just seconds per period.  The
    % inverse is the sinewave frequency in period per seconds i.e. the "hertz"
    % cycles per second.

end

