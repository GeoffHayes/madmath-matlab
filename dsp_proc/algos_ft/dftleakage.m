%**************************************************************************
% Plots the DFT leakage given the different input sinusoid cycle(s), an
% N-point DFT, and the sampling frequency.  It is assumed that each
% sinusoid has unity amplitude.
%
% @param   k   The input sinusoid cycles.
% @param   fs  The sampling frequency (Hz).
% @param   N   The size of the DFT.
% @param   mag The magnitude output of the DFT.
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
function dftleakage(k,fs,N,mag)

    N = abs(N(1));
    
    % frequency resolution
    fres = fs/N;
    
    numRowPlots = length(k);

    figure;
    
    % we only plot the first N/2 bins as all inputs to the DFT are ASSUMED
    % TO BE real and so the first N/2 bins are independent
    
    for i=1:numRowPlots

        % determine the frequencies given the frequency resolution and the
        % number of bins
        f = ((0:N/2)*fres)';
        
        % apply the sinc function
        sincf = (N/2)*sinc(k(i)-f);
        
        % construct an interval of data for the continuous positive
        % frequency spectrum of a discrete cosine sequence
        t=(0:0.1:N/2)';
        
        % apply the sinc function 
        sinct = (N/2)*sinc(k(i)/fres-t);

        subplot(numRowPlots,2,2*(i-1)+1);
        hold on;
        plot(f/fres,sincf,'rs');
        plot(t,sinct,'b');
        line([t(1) t(end)],[0 0],'Color','k');
        
        title(['Amplitude response as a function of bin index m (k=' ...
            num2str(k(i)/fres) ')']);
        xlabel('Frequency (m)');
        axis([0 max(t) min(min(sincf),min(sinct)) max(max(sincf),max(sinct))]);
        
        subplot(numRowPlots,2,2*i);
        hold on;
        plot((0:N/2),mag(1:N/2+1),'rs');
        plot(t,abs(sinct),'b');
        line([t(1) t(end)],[0 0],'Color','k');
        
        title('Magnitude response as a function of frequency (Hz)');
        xlabel('Frequency (kHz)');   
        axis([0 max(t) min(min(sincf),min(sinct)) max(max(sincf),max(sinct))]);
 
    end
   
end

