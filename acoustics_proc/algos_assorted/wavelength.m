%**************************************************************************
% Calculates the wavelength given the frequency and speed of sound.
%
% @param   f      The frequency of the sound (Hz).
% @param   c      The speed of sound (m/s).
%
% @return  The wavelength.
%
% [1] http://en.wikipedia.org/wiki/Wavelength
%**************************************************************************
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
function [lambda] = wavelength(f,c)

    if ~exist('c','var')
        % use default speed of sound is 343 m/s (at room temperature and  
        % atmospheric pressure)
        c = 343;
    end
    
    lambda = c/f;
    
end

