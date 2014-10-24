%**************************************************************************
% Rectangular window function (also known as the uniform or boxcar window):
% 
%               w(n) = 1, for n=0,1,2,....,N-1
%
% @param   N   The length of the window.
%
% @return  The rectangular window vector.
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
function [w] = wdw_rect(N)

    % Windows are designed to reduce the problem of spectral leakage. 
    % Leakage occurs with every frequency component that is not harmonic 
    % with the FFT fundamental. The correlation for such an inharmonic 
    % sinusoid will spread over all bins. There is a peak and a base in the 
    % spectrum of this sinusoid, but the contour can be manipulated by 
    % applying one or the other window type. Window functions are 
    % multiplied with the FFT input signal. They modulate that signal, 
    % thereby altering the leakage character of the analysis. 
    % http://www.katjaas.nl/FFTwindow/FFTwindow.html
    
    % The rectangular window has excellent resolution characteristics for 
    % sinusoids of comparable strength, but it is a poor choice for 
    % sinusoids of disparate amplitudes. This characteristic is sometimes 
    % described as low-dynamic-range.
    
    N = abs(floor(N(1)));
    
    w=ones(N,1);
    
    % the frequency response of the retangular window, w(m), is as follows, 
    % where F is the Fourier Transform and int(a,b) is the integral of the
    % function integrated with respect to x within the bounds a and b:
    %
    %   F{w(x)} = int(-b,b)exp(-i*2*pi*x*m)dx
    %           = (exp(-i*2*pi*m*b) - exp(i*2*pi*m*b))/(-i*2*pi*m)
    %           = (cos(2*pi*m*b) - i*sin(2*pi*m*b) - cos(2*pi*m*b) -
    %                   i*sin(2*pi*m*b))/(-i*2*pi*m)
    %           = -i(2*sin(2*pi*m*b)/(-i*2*pi*m)
    %           = sin(2*pi*m*b)/(pi*m)
    %           = F(m)
    %
    % if the interval that we are integrating over is of length N (the FFT
    % size) then b=N/2 which gives:
    %
    %   F(m) = F{w(x)} = sin(2*pi*m*N/2)/(pi*m) 
    %                  = sin(pi*m*N)/(pi*m)
    
    % to plot the frequency response of the rectangular window:
    % x2 = @(t)(sin(1024*pi*t)./(pi*t));
    % close all;t=0:0.001:1;W=x2(t);W(1)=ceil(W(2));plot(t,todB(W,1))
    
    % note that F(w(m)) = sin(pi*m*N)/(pi*m) =
    % N/2*sin(pi*m)*cos(pi*m)/(pi*m) :=: N/2*sin(pi*m)/(pi*m) which leads
    % to the alternative:
    %
    % x2 = @(t)(1024/2*sinc(t));
    % close all;t=-200:0.001:200;W=x2(t);W(W==0)=max(W)+eps;plot(t,todB(W,1))
    
%     code to simulate the fourier transform of the rectangular window
%     N = 1024;
%     res = 0.001;
%     F   = zeros(1/res-1,1);
%     at  = 1;
%     for m=-1:res:1
%         % ignore the DC component
%         if m~=0
%             func  = @(x)(exp(-1i*2*pi*x*m));
%             F(at) = integral(func,-N/2,N/2);
%             at    = at + 1;
%         end
%     end

end