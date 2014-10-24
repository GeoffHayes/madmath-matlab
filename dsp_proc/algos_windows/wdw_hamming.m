%**************************************************************************
% Hamming window function:
% 
%               w(n) = 0.54-0.46cos(2pin/(N-1)), n=0,1,....,N-1
%
% @param   N   The length of the window.
%
% @return  The Hamming window vector.
%
% http://www.mathworks.com/help/signal/ref/hamming.html
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
function [w] = wdw_hamming(N)

    % Windows are designed to reduce the problem of spectral leakage. 
    % Leakage occurs with every frequency component that is not harmonic 
    % with the FFT fundamental. The correlation for such an inharmonic 
    % sinusoid will spread over all bins. There is a peak and a base in the 
    % spectrum of this sinusoid, but the contour can be manipulated by 
    % applying one or the other window type. Window functions are 
    % multiplied with the FFT input signal. They modulate that signal, 
    % thereby altering the leakage character of the analysis. 
    % http://www.katjaas.nl/FFTwindow/FFTwindow.html

    N = abs(floor(N(1)));
        
    w=zeros(N,1);
    
    if N==1
        w = 1;
    else
        for n=0:N-1        
            w(n+1) = 0.54 - 0.46*cos(2*pi*n/(N-1));
        end
    end

end