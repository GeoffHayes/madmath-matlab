%**************************************************************************
% Performs a transformation on the independent variable n of the
% discretized data:
%
%   1.  n-n0 : y[n] = x[n-n0] indicates that the output is shifted by n0
%       samples to the right (this is a delay).
%   2.  n+n0 : y[n] = x[n+n0] indicates that the output is shifted by n0
%       samples to the left (this is an advance).
%   3.  -n   : y[n] = x[-n] indicates that the output is reversed/flipped.
%   4.  Mn   : y[n] = x[Mn] with M a positive integer indicates that every
%       Mth sample is taken of x[n] (down-sampling).
%   5.  n/M  : y[n] = x[n/M] with M a positive integer indicates that y[n]
%       = x[n/N] for n=...,-2N,-N,0,N,2N,... and 0 otherwise (up-sampling).
%
% @param   plotSec       The period in seconds to plot the data (can be 
%                        interval).
% @param   fs            The sampling rate in samples per second.
% @param   wavFnc        The wave function.
% @param   tfrmFncs      The transform function(s).
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
function wavtfrma2d(plotSec,fs,wavFnc,tfrmFncs)

    % use a gap of 1/1,000th of a second
    analogGap = 0.00001;

    % if the sampling rate is too great, fix to the max allowed (which is the
    % gap)
    numSamplesPerSecond = 1/analogGap;
    if fs > numSamplesPerSecond
        fs = numSamplesPerSecond;
    end    
    ts = 1/fs;
    
    if length(plotSec)==1
        T = (0:analogGap:plotSec)';
        N = (0:ts:plotSec)';
    else
        if plotSec(1)>plotSec(2)
            temp = plotSec(2);
            plotSec(2) = plotSec(1);
            plotSec(1) = temp;
        end
        T = (plotSec(1):analogGap:plotSec(2))';
        N = (plotSec(1):ts:plotSec(2))';
    end
    
    % given the time, compute the input analog time domain data evaluated
    % for the input function
    anaTimeDomainData = wavFnc(T);
    
    % apply the A/D converter against the analog time domain data
    digTimeDomainData = wavFnc(N);
    
    % apply the transforms
    numTfrms               = length(tfrmFncs);
    digTimeDomainTfrmdData = zeros(length(N),numTfrms);
    Np                     = zeros(length(N),numTfrms);
    
    for i=1:numTfrms
       Np(:,i) = tfrmFncs{i}(N);
       digTimeDomainTfrmdData(:,i) = wavFnc(Np(:,i));
    end

    maxYAxis = max([max(anaTimeDomainData)  max(digTimeDomainData) ...
                    max(digTimeDomainTfrmdData(:))]);
                
    maxXAxis = max([max(T)  max(N) ...
                    max(Np(:))]);                
    
    figHand = figure;
    
    set(figHand,'Name','Contrasting Discrete Time Domain Transforms');
    
    numSubPlots = 2 + numTfrms;  % 2 for the analog and digital equivalent


    % plot the waves individually and then the combined wave
    subplot(numSubPlots,1,1);
    plot(T,anaTimeDomainData(:,1));
    xlabel('Time (sec)');
    line([T(1) T(end)],[0 0],'Color','k');
    tlt1a = char(wavFnc);
    tlt1a = tlt1a(5:end);
    title(['x_1(t) = ' tlt1a]);
    axis([-maxXAxis maxXAxis -maxYAxis maxYAxis]);
    anaXTicks = get(gca,'XTickLabel');

    % draw a lollipop plot for the Digital data
    subplot(numSubPlots,1,2);
    plot(N,digTimeDomainData(:,1),'r.', 'MarkerSize',8);
    xlabel('Time (sec)');
    line([T(1) T(end)],[0 0],'Color','k');
    title('y_1[n]');
    axis([-maxXAxis maxXAxis -maxYAxis maxYAxis]);   
    for i=1:length(N)
       line([N(i) N(i)],[0 digTimeDomainData(i,1)],'Color','r'); 
    end
    
    for i=1:numTfrms

        % draw a lollipop plot for the Digital data
        subplot(numSubPlots,1,2+i);
        plot(N,digTimeDomainTfrmdData(:,i),'r.', 'MarkerSize',8);
        xlabel('Time (sec)');
        line([T(1) T(end)],[0 0],'Color','k');
        title(['y_' num2str(i+1) '[n]']);
        axis([-maxXAxis maxXAxis -maxYAxis maxYAxis]);   
        for j=1:length(N)
           line([N(j) N(j)],[0 digTimeDomainTfrmdData(j,i)],'Color','r'); 
        end
    end

end


% usage: creates a single triangle with height 3 and centre/peak at 3 and
% then performs the various shifts against the data
%
% x1=@(t)singletri(3,3,1,t);
% t1=@(n)tfrm_shift(n,-2,x1)
% t2=@(n)tfrm_flip(n,x1)
% t3=@(n)tfrm_downsample(n,2,x1)
% t4=@(n)tfrm_upsample(n,2,x1)
% close all;wavtfrma2d([-10 10],1,x1,{t1 t2 t3 t4})
