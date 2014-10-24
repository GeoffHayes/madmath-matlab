%**************************************************************************
% Power of cosine window function:
% 
%               w(n) = cos^a((pi*n)/(N-1) - pi/2)
%
% @param   N   The length of the window.
% @param   a   The cosine power (positive integer).
%
% @return  The power of cosine window vector.
%
% [1] http://en.wikipedia.org/wiki/Window_function
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
function [w] = wdw_cos(N,a)

    if ~exist('a','var')
        a = 1;
    else
        a = abs(floor(a(1)));
    end
    
    N = abs(floor(N(1)));

    % Windows are designed to reduce the problem of spectral leakage. 
    % Leakage occurs with every frequency component that is not harmonic 
    % with the FFT fundamental. The correlation for such an inharmonic 
    % sinusoid will spread over all bins. There is a peak and a base in the 
    % spectrum of this sinusoid, but the contour can be manipulated by 
    % applying one or the other window type. Window functions are 
    % multiplied with the FFT input signal. They modulate that signal, 
    % thereby altering the leakage character of the analysis. 
    % http://www.katjaas.nl/FFTwindow/FFTwindow.html
    
    w=zeros(abs(N),1);
    
    for n=0:N-1
       w(n+1) = (cos((pi*n )/(N-1) - pi/2))^a;        
    end
    
    

end