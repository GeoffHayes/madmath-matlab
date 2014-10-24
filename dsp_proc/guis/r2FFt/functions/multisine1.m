function [y,Fs] = multisine1(t)

Fs  = length(t);

y = 100*sin(2*pi*50*t)+1*sin(2*pi*123*t)+2*sin(2*pi*203*t)+ ...
    3*sin(2*pi*223*t)+4*sin(2*pi*331*t)+5*sin(2*pi*2812*t)+ ...
    6*sin(2*pi*5752*t)+7*sin(2*pi*7993*t);

% see http://techteach.no/simview/aliasing/index.php for aliasing

% Sampling av continuous-time signals takes place in all applications where 
% a computer is used to read measurement data. If the sampling rate is too 
% low compared to the frequency of the signal to be sampled, aliasing
% occurs. Aliasing means that the discrete-time signal gets a lower 
% frequency than the original signal. Obviously, sampling can cause 
% problems in applications, e.g. in audio applications.

% Given a continuous-time signal, x(t), of frequency fcont [Hz] which is 
% sampled with sampling frequency fs [Hz]. The Nyquist frequency fN is 
% defined as half of the sampling frequency:

% fN = fs/2

% It ca be shown that if the signal frequency fcont is larger than fN, 
% there is aliasing, which implies that the resulting discrete-time signal, 
% xs, gets a frequency, fdisc, which is smaller than the signal frequency 
% fcont. Figure 1 shows the relation between fdisc and fcont.