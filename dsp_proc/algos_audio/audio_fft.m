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
function audio_fft(filename,fftsize)

    % read the audio data
    try
        [Z,fs] = audioread(filename);
    catch
        err('audio_fft: file %s cannot be found',filename);
    end
    
    totalSamples = length(Z);
    atSample     = 1;
    N            = 2^nextpow2(fftsize);
    
    while atSample < totalSamples
        
        % grab a block of data
        if atSample+N>totalSamples
            x = Z(atSample:end);
        else
            x = Z(atSample:atSample+N);
        end
        atSample = atSample + N;
        
        % do the FFT
        X = r2fft(x,N);
  
    end

end

