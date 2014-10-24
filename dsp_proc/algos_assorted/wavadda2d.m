%**************************************************************************
% Adds/combines (any) two waves. The input wave (x[n]) is plotted,  
% followed by the output wave (y(x[n]), followed by the discretized 
% digital equivalent of the latter (outptu) data, simulating the 
% analog-to-digital (A/D) conversion.
%
% @param   plotSec       The period in seconds to plot the data.
% @param   fs            The sampling rate in samples per second.
% @param   inputFncsX    The input functions for each wave.
% @param   outputFncsY   The output functions (which may be linear or
%                        non-linear) for each wave.
%
% @throw   Error if the input vectors are not of the same length.
% @throw   Error if the input vectors are not of size two.
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
function wavadda2d(plotSec,fs,inputFncsX,outputFncsY)

    MAX_NUM_WAVES=2;
    
    % use a gap of 1/1,000th of a second
    analogGap = 0.00001;

    time = (0:analogGap:plotSec)';
    
    if length(inputFncsX)~=length(outputFncsY)
        error('wavedda2d: input vectors are of different lengths');
    elseif length(inputFncsX)~=MAX_NUM_WAVES
        error('wavedda2d: input vectors are not of size two');
    end
    
    % if the sampling rate is too great, fix to the max allowed (which is the
    % gap)
    numSamplesPerSecond = 1/analogGap;
    if fs > numSamplesPerSecond
        fs = numSamplesPerSecond;
    end

    ts = 1/fs;
    
    % given the time, compute the input analog time domain data evaluated
    % for all functions X
    anaTimeDomainDataX = zeros(length(time),MAX_NUM_WAVES);
    
    for i=1:length(inputFncsX)
        anaTimeDomainDataX(:,i) = inputFncsX{i}(time);
    end
            
    % given the input analog time domain data, compute the output analog 
    % time domain data evaluated for all functions Y
    anaTimeDomainDataY = zeros(length(time),MAX_NUM_WAVES);
    
    for i=1:length(outputFncsY)
        anaTimeDomainDataY(:,i) = outputFncsY{i}(anaTimeDomainDataX(:,i));
    end
       
    % apply the A/D converter against the analog time domain data
    N                  = (0:ts:plotSec)';  % labelled as N for y[n], where n belongs to N
    digTimeDomainDataX = zeros(length(N),MAX_NUM_WAVES);
    
    for i=1:length(inputFncsX)
        digTimeDomainDataX(:,i) = inputFncsX{i}(N);
    end
    
    digTimeDomainDataY = zeros(length(N),MAX_NUM_WAVES);
    
    for i=1:length(outputFncsY)
        digTimeDomainDataY(:,i) = outputFncsY{i}(digTimeDomainDataX(:,i));
    end
    
    % now combine the analog output data (X)
    anaCombndTimeDomainDataX = anaTimeDomainDataX(:,1) + ...
        anaTimeDomainDataX(:,2);    
    
    % now combine the analog output data (Y) using the first output
    % function Y
    anaCombndTimeDomainDataY = outputFncsY{1}(anaTimeDomainDataX(:,1) + ...
        anaTimeDomainDataX(:,2));
    
    % now combine the digital output data (Y) using the first output
    % function Y
    digCombndTimeDomainDataY = outputFncsY{1}(digTimeDomainDataX(:,1) + ...
        digTimeDomainDataX(:,2));
    
    
    maxYAxis = max([max(anaTimeDomainDataX)  max(anaTimeDomainDataY) ...
                    max(anaCombndTimeDomainDataY) max(digCombndTimeDomainDataY)]);
    
    figHand = figure;
    
    set(figHand,'Name',['Contrasting Time Domain Waves (x[n] vs' ...
            ' y[n]) vs A/D of y[n])']);


    % plot the waves individually and then the combined wave
    subplot(3,3,1);
    plot(time,anaTimeDomainDataX(:,1));
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    tlt1a = char(inputFncsX{1});
    tlt1a = tlt1a(5:end);
    title(['x_1(t) = ' tlt1a]);
    axis([0 time(end) -maxYAxis maxYAxis]);
    anaXTicks = get(gca,'XTickLabel');
    
    subplot(3,3,2);
    plot(time,anaTimeDomainDataY(:,1));
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    tlt1b = char(outputFncsY{1});
    tlt1b = strrep(tlt1b(5:end),'y','(x_1(t))');
    title(['y_1(t) = ' tlt1b]);
    axis([0 time(end) -maxYAxis maxYAxis]);
    
    % draw a lollipop plot for the Digital data
    subplot(3,3,3);
    plot(N,digTimeDomainDataY(:,1),'r.', 'MarkerSize',8);
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    title('y_1[n]');
    axis([0 time(end) -maxYAxis maxYAxis]);   
    for i=1:length(N)
       line([N(i) N(i)],[0 digTimeDomainDataY(i,1)],'Color','r'); 
    end
    set(gca,'XTickLabel',anaXTicks);

    subplot(3,3,4);
    plot(time,anaTimeDomainDataX(:,2));
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    tlt2a = char(inputFncsX{2});
    tlt2a = tlt2a(5:end);
    title(['x_2(t) = ' tlt2a]);      
    axis([0 time(end) -maxYAxis maxYAxis]);
    
    subplot(3,3,5);
    plot(time,anaTimeDomainDataY(:,2));
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    tlt2b = char(outputFncsY{2});
    tlt2b = strrep(tlt2b(5:end),'y','(x_2(t))');
    title(['y_2(t) = ' tlt2b]);
    axis([0 time(end) -maxYAxis maxYAxis]);
    
    % draw a lollipop plot for the Digital data
    subplot(3,3,6);
    plot(N,digTimeDomainDataY(:,2),'r.', 'MarkerSize',8);
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    title('y_2[n]');
    axis([0 time(end) -maxYAxis maxYAxis]);   
    for i=1:length(N)
       line([N(i) N(i)],[0 digTimeDomainDataY(i,2)],'Color','r'); 
    end
    set(gca,'XTickLabel',anaXTicks);    
    
    subplot(3,3,7);
    plot(time,anaCombndTimeDomainDataX);
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    title('x_3(t) = x_1(t) + x_2(t)');    
    axis([0 time(end) -maxYAxis maxYAxis]);
    
    subplot(3,3,8);
    plot(time,anaCombndTimeDomainDataY);
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    temp = char(outputFncsY{1});
    temp = strrep(temp(5:end),'y','(x_1(t)+x_1(t))');
    title(['y_3(t)=' temp]);
    axis([0 time(end) -maxYAxis maxYAxis]);
    
    % draw a lollipop plot for the Digital data
    subplot(3,3,9);
    plot(N,digCombndTimeDomainDataY,'r.', 'MarkerSize',8);
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    title('y_3[n]');
    axis([0 time(end) -maxYAxis maxYAxis]);   
    for i=1:length(N)
       line([N(i) N(i)],[0 digCombndTimeDomainDataY(i)],'Color','r'); 
    end
    set(gca,'XTickLabel',anaXTicks);       
    
end

% usage:
%
% x1=@(t)1.0*sin(2*pi*1*t)
% x2=@(t)1.0*sin(2*pi*3*t)
% y1=@(y)y.^2
% y2=@(y)y.^2
% 
% close all;wavadda2d(1,25,{x1 x2},{y1 y2})
