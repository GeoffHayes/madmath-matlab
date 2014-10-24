%**************************************************************************
% Adds/combines two sine waves, plotting the original sine wave for the
% given frequency; the output of the linear system which is a function of
% the input sine wave, the frequency scalar, and amplitude; and the
% combined sine waves.
%
% @param   plotSec       The period in seconds to plot the data.
% @param   frequenciesHz The frequency of each sine wave in hertz (number
%                        of cycles of the sine wave per second).
% @param   freqScalars   The frequency scalars for each sine wave.
% @param   amplitudes    The amplitudes for each sine wave.
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
function sinadd(plotSec,frequenciesHz,freqScalars,amplitudes)

    % use a gap of 1/1,000th of a second
    gap = 0.001;

    time = (0:gap:plotSec)';
    
    if length(freqScalars)~=length(amplitudes) && ...
            length(frequenciesHz)~=length(amplitudes)
        error('sinadd: input vectors are of different lengths');
    elseif length(freqScalars)~=2
        error('sinadd: input vectors are not of size two');
    end

    % calculate the analog data using the time-domain equation with the
    % application of the freqency and amplitude scalars
    scaledFrequencies = frequenciesHz.*freqScalars;
    timeDomainDataB = sin(2*pi*time*scaledFrequencies);
    for i=1:length(amplitudes)
        timeDomainDataB(:,i) = timeDomainDataB(:,i)*amplitudes(i);
    end
    
    % calculate the analog data using the time-domain equation with only
    % difference being the frequency
    timeDomainDataA = sin(2*pi*time*scaledFrequencies);    
        
    
    cTimeDomainDataA = timeDomainDataA(:,1) + timeDomainDataA(:,2);
    cTimeDomainDataB = timeDomainDataB(:,1) + timeDomainDataB(:,2);
    
    maxYAxis = max([max(timeDomainDataA)  max(timeDomainDataB) ...
                    max(cTimeDomainDataA) max(cTimeDomainDataB)]);
    
    figHand = figure;
    set(figHand,'Name',['Contrasting Time Domain Waveforms (Original vs' ...
        ' Linear System Output)']);

    % plot the sine waves individually and then the combined wave
    subplot(3,2,1);
    plot(time,timeDomainDataA(:,1));
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    tlt1a = 'sin(2pif1nt)';
    title([tlt1a ' with f1=' num2str(frequenciesHz(1)) ' Hz']);
    axis([0 time(end) -maxYAxis maxYAxis]);
    
    subplot(3,2,2);
    plot(time,timeDomainDataB(:,1));
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    tlt1b = [num2str(amplitudes(1)) 'sin(2pi' num2str(freqScalars(1)) ...
        'f1nt)'];
    title([tlt1b ' with f1=' num2str(frequenciesHz(1)) ' Hz']);
    axis([0 time(end) -maxYAxis maxYAxis]);

    subplot(3,2,3);
    plot(time,timeDomainDataA(:,2));
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    tlt2a = 'sin(2pif1nt)';
    title([tlt2a ' with f1=' num2str(frequenciesHz(2)) ' Hz']);   
    axis([0 time(end) -maxYAxis maxYAxis]);
    
    subplot(3,2,4);
    plot(time,timeDomainDataB(:,2));
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    tlt2b = [num2str(amplitudes(2)) 'sin(2pi' num2str(freqScalars(2)) ...
        'f2nt)'];
    title([tlt2b ' with f2=' num2str(frequenciesHz(2)) ' Hz']);
    axis([0 time(end) -maxYAxis maxYAxis]);
    
    subplot(3,2,5);
    plot(time,cTimeDomainDataA);
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    title([tlt1a ' + ' tlt2a]);    
    axis([0 time(end) -maxYAxis maxYAxis]);
    
    subplot(3,2,6);
    plot(time,cTimeDomainDataB);
    xlabel('Time (sec)');
    line([time(1) time(end)],[0 0],'Color','k');
    title([tlt1b ' + ' tlt2b]);
    axis([0 time(end) -maxYAxis maxYAxis]);
    
    % what we are trying to show is that there are two frequency components
    % in the combined equation (waveform) that of frequency1 and 
    % frequency2.  In the frequency domain, we should be able to pull this 
    % frequency information out of the time-domain equation.
    
end

% sinadd(5,[1 3],[1 1],[-0.5 -0.5])
