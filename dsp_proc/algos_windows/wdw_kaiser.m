%**************************************************************************
% Kaiser window function:
% 
%               w(n) = Io{pi*a*sqrt(1-((2n)/(N-1) - 1)^2}/Io{pi*a}, 
%               for n=0,1,....,N-1
%
% @param   N   The length of the window.
% @param   a   The non-negative real number that determines the shape of
%              the window.  In the frequency domain, it determines the
%              trade-off between main lobe width and side lobe level.
%
% @return  The Kaiser window vector.
%
% http://en.wikipedia.org/wiki/Kaiser_window
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
function [w] = wdw_kaiser(N,a)

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
    a = abs(a(1));
    
    w = zeros(N,1);
    
    for i=0:N-1
        w(i+1) = besseli(0,pi*a*sqrt(1-((2*i/(N-1))-1)^2))/...
            besseli(0,pi*a);
    end

end