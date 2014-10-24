function [y,Fs] = guitar1sec(~)

persistent pFs;
persistent data;
persistent atSample;

if isempty(data)
    [data,pFs] = audioread('geoffEstring.wav');
    atSample  = 1;
end

numSamples = length(data);
strtSample = atSample;
endSample  = atSample+pFs-1;

if endSample>numSamples
    y = zeros(pFs,1);
    endSample = min(endSample,length(data));
    y(1:endSample-strtSample+1) = data(strtSample:end);
    atSample = 1;
else
    y        = data(atSample:atSample+pFs-1,1);
    atSample = atSample + pFs;
end

Fs       = pFs;

sound(y,Fs);

