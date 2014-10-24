%**************************************************************************
% Adds/combines two sine waves, plotting the original sine wave for the
% given frequency; the output of the linear system which is a function of
% the input sine wave, the frequency scalar, and amplitude; and the
% combined sine waves.  The discretized (A/D of the LSO) is plotted for
% each input and (LSO) output pair.  Note that if a non-linear function
% is provided, then the Non-Linear System Output (NLSO) is plotted for the 
% analog and the digital waveforms.
%
% @param   plotSec       The period in seconds to plot the data.
% @param   fs            The sampling rate in samples per second.
% @param   frequenciesHz The frequency of each sine wave in hertz (number
%                        of cycles of the sine wave per second).
% @param   freqScalars   The frequency scalars for each sine wave.
% @param   amplitudes    The amplitudes for each sine wave.
% @param   nonLinFnc     A non-linear function that can be applied to the
%                        system.
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
function sinadda2d(plotSec,fs,frequenciesHz,freqScalars,amplitudes,nonLinFnc)

    % use a gap of 1/1,000th of a second
    analogGap = 0.00001;

    time = (0:analogGap:plotSec)';
    
    if length(freqScalars)~=length(amplitudes) && ...
            length(frequenciesHz)~=length(amplitudes)
        error('sinadda2d: input vectors are of different lengths');
    elseif length(freqScalars)~=2
        error('sinadda2d: input vectors are not of size two');
    end
    
    if ~exist('nonLinFnc','var')
        nonLinFnc = {};
    end
    
    % if the sampling rate is too great, fix to the max allowed (which is the
    % gap)
    numSamplesPerSecond = 1/analogGap;
    if fs > numSamplesPerSecond
        fs = numSamplesPerSecond;
    end

    ts = 1/fs;

    % calculate the analog data using the time-domain equation with the
    % application of the freqency and amplitude scalars
    scaledFrequencies = frequenciesHz.*freqScalars;
    timeDomainDataB = sin(2*pi*time*scaledFrequencies);
    
    % apply the A/D converter against the sine wave (analog data)
    digTime            = (0:ts:plotSec)';                 
    digTimeDomainDataB = sin(2*pi*digTime*scaledFrequencies);       
    
    for i=1:length(amplitudes)
        timeDomainDataB(:,i)    = timeDomainDataB(:,i)*amplitudes(i);
        digTimeDomainDataB(:,i) = digTimeDomainDataB(:,i)*amplitudes(i);    
    end
    
    % calculate the analog data using the time-domain equation with the
    % application of the freqency scalar
    timeDomainDataA = sin(2*pi*time*scaledFrequencies);    

    cTimeDomainDataA    = timeDomainDataA(:,1) + timeDomainDataA(:,2);

    if isempty(nonLinFnc)
        cTimeDomainDataB    = timeDomainDataB(:,1) + ...
            timeDomainDataB(:,2);
        cDigTimeDomainDataB = digTimeDomainDataB(:,1) + ...
            digTimeDomainDataB(:,2);
    else
        cTimeDomainDataB    = nonLinFnc{1}(timeDomainDataB(:,1) + ...
            timeDomainDataB(:,2));
        cDigTimeDomainDataB = nonLinFnc{1}(digTimeDomainDataB(:,1) + ...
            digTimeDomainDataB(:,2));
        
        for i=1:length(amplitudes)
            timeDomainDataB(:,i)    = nonLinFnc{1}(...
                timeDomainDataB(:,i)*amplitudes(i));
            digTimeDomainDataB(:,i) = nonLinFnc{1}(...
                digTimeDomainDataB(:,i)*amplitudes(i));    
        end
    end
    
    maxYAxis = max([max(timeDomainDataA)  max(timeDomainDataB) ...
                    max(cTimeDomainDataA) max(cTimeDomainDataB)]);
    
    figHand = figure;
    
    if isempty(nonLinFnc)
        set(figHand,'Name',['Contrasting Time Domain Waveforms (Original vs' ...
            ' Linear System Output (LSO) vs A/D of LSO)']);
    else
        set(figHand,'Name',['Contrasting Time Domain Waveforms (Original vs' ...
        ' Non-Linear System Output (NLSO) vs A/D of NLSO)']);
    end

    % plot the sine waves individually and then the combined wave
    subplot(3,3,1);
    plot(time,timeDomainDataA(:,1));
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    tlt1a = 'sin(2pif1nt)';
    title([tlt1a ' with f1=' num2str(frequenciesHz(1)) ' Hz (A)']);
    axis([0 time(end) -maxYAxis maxYAxis]);
    anaXTicks = get(gca,'XTickLabel');
    
    subplot(3,3,2);
    plot(time,timeDomainDataB(:,1));
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    tlt1b = [num2str(amplitudes(1)) 'sin(2pi' num2str(freqScalars(1)) ...
        'f1nt)'];
    title([tlt1b ' with f1=' num2str(frequenciesHz(1)) ' Hz (A)']);
    axis([0 time(end) -maxYAxis maxYAxis]);
    
    % draw a lollipop plot for the Digital data
    subplot(3,3,3);
    plot(digTime,digTimeDomainDataB(:,1),'r.', 'MarkerSize',8);
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    title([tlt1b ' with f1=' num2str(frequenciesHz(1)) ' Hz (D)']);
    axis([0 time(end) -maxYAxis maxYAxis]);   
    for i=1:length(digTime)
       line([digTime(i) digTime(i)],[0 digTimeDomainDataB(i,1)],'Color','r'); 
    end
    set(gca,'XTickLabel',anaXTicks);

    subplot(3,3,4);
    plot(time,timeDomainDataA(:,2));
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    tlt2a = 'sin(2pif1nt)';
    title([tlt2a ' with f1=' num2str(frequenciesHz(2)) ' Hz (A)']);   
    axis([0 time(end) -maxYAxis maxYAxis]);
    
    subplot(3,3,5);
    plot(time,timeDomainDataB(:,2));
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    tlt2b = [num2str(amplitudes(2)) 'sin(2pi' num2str(freqScalars(2)) ...
        'f2nt)'];
    title([tlt2b ' with f2=' num2str(frequenciesHz(2)) ' Hz (A)']);
    axis([0 time(end) -maxYAxis maxYAxis]);
    
    % draw a lollipop plot for the Digital data
    subplot(3,3,6);
    plot(digTime,digTimeDomainDataB(:,2),'r.', 'MarkerSize',8);
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    title([tlt1b ' with f1=' num2str(frequenciesHz(2)) ' Hz (D)']);
    axis([0 time(end) -maxYAxis maxYAxis]);   
    for i=1:length(digTime)
       line([digTime(i) digTime(i)],[0 digTimeDomainDataB(i,2)],'Color','r'); 
    end
    set(gca,'XTickLabel',anaXTicks);    
    
    subplot(3,3,7);
    plot(time,cTimeDomainDataA);
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    title([tlt1a ' + ' tlt2a ' (A)']);    
    axis([0 time(end) -maxYAxis maxYAxis]);
    
    subplot(3,3,8);
    plot(time,cTimeDomainDataB);
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    title([tlt1b ' + ' tlt2b ' (A)']);
    axis([0 time(end) -maxYAxis maxYAxis]);
    
    % draw a lollipop plot for the Digital data
    subplot(3,3,9);
    plot(digTime,cDigTimeDomainDataB,'r.', 'MarkerSize',8);
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    title([tlt1b ' + ' tlt2b ' (D)']);
    axis([0 time(end) -maxYAxis maxYAxis]);   
    for i=1:length(digTime)
       line([digTime(i) digTime(i)],[0 cDigTimeDomainDataB(i)],'Color','r'); 
    end
    set(gca,'XTickLabel',anaXTicks);       
    
end

% sinadda2d(1,40,[1 3],[1 1],[-0.5 1.5])
% sinadda2d(1,10,[1 3],[1 1],[-0.5 -0.5])

