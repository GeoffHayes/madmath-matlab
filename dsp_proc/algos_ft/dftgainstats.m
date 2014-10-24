%**************************************************************************
% Generates some plots illustrating the effects of the number of bins (FFT
% size) against a signal hidden in noise.  Here we are using the DFT to
% detect signal energy embedded, and this is referred to as the processing
% gain because the DFT can pull signals out of background noise ([1]).
%
% @param   func   The function used to generate a signal embedded in noise.
% @param   fs     The sampling frequency (samples per second) that the data
%                 will be generated at.
% @param   N      One or more FFT block sizes.  As the block size 
%                 increases, the signal hidden in the noise is more clear.
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
function dftgainstats(func,fs,N)

    % we sort from smallest N(i) to largest N(i) to demonstrate that the
    % DFT's output SNR increases as N gets larger because a DFT's output
    % noise standard deviation (rms) value is proportional to aqrt(N), and
    % the DFT's output magnitude for the bin containing the signal tone is
    % proportional to N.  According to [1], for real inputs, if N>Np, an
    % N-point DFT's output SNR_N increases over the Np-point DFT SNR_Np by
    % the following relationship:
    %
    % SNR_N=SNR_Np + 20*log10(sqrt(N/Np))
    %
    % By the above, if N=2Np, then the DFT's output SNR increases by 3 dB.
    N = sort(ceil(abs(N)));
    
    numPlots = length(N);
    
    data = {};
    
    for i=1:numPlots
        
        % generate the analog and discrete data
        [~,data{i}] = wava2d(N(i),fs,func); %#ok<AGROW>
        
    end
    
    figure;    
    
    for i=1:numPlots
        
        % perform the DFT given the block size
        X     = dft(data{i},N(i));     
        
        % convert the frequency domain data to power
        P = abs(X).^2;
        
        % plot the power data
        subplot(numPlots,1,i);
        
        plot(todB(P(1:N(i)/2),1),'b');
        title(['DFT Processing Gain with N=' num2str(N(i))]);
        xlabel('DFT bin number');
        ylabel('Bin power in dB');
        
    end
end

% usage:
% x1 = @(t)(sin(2*pi*20*t) + 2.0*randn(size(t)));
% dftgainstats(x1,64,[64 256 1024])

