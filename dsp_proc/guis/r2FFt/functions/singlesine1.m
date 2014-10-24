function [y,Fs] = singlesine1(t)

amp = 12.0;
frq = 75.0;
Fs  = length(t);

y = amp*sin(2*pi*frq*t);

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

% It can be shown that if the signal frequency fcont is larger than fN, 
% there is aliasing, which implies that the resulting discrete-time signal, 
% xs, gets a frequency, fdisc, which is smaller than the signal frequency 
% fcont. Figure 1 shows the relation between fdisc and fcont.