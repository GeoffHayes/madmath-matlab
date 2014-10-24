%**************************************************************************
% Triangular window function:
% 
%               w(n) = n/(N/2),   for n=0,1,2,....,N/2
%                    = 2-n/(N/2), for n=N/2+1,....,N-1
%
% @param   N   The length of the window.
%
% @return  The triangular window vector.
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
function [w] = wdw_triangle(N,x)

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

    if mod(N,2)~=0
        N = N + 1;
    end
    
    if ~exist('x','var')

        w=zeros(abs(N),1);

        for n=0:N/2        
            w(n+1) = n/(N/2);
        end

        for n=N/2+1:N-1
            w(n+1) = (N-n)/(N/2);
        end
        
    else
        w = zeros(size(x));
        for i=1:length(w)
            if x(i)<=0
                w(i) = 1+x(i)/(N/2);
            elseif x(i)<=(N/2-1)
                w(i) = 1-x(i)/(N/2);
            end
        end
    end
    
    
%     code to simulate the fourier transform of the triangular window (not
%     really)
%     N = 1024;
%     res = 0.01;
%     F   = zeros(1/res-1,1);
%     at  = 1;
%     for m=-1:res:1
%         % ignore the DC component
%         %if m~=0
%             func  = @(x)(wdw_triangle(N,x).*exp(-1i*2*pi*x*m));
%             F(at) = integral(func,-N/2,N/2);
%             at    = at + 1;
%         %end
%     end

end