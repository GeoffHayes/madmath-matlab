%**************************************************************************
% Converts the (DFT) frequency domain sequence power spectrum sequence in 
% decibels.
%
%           P(m) = 10*log10(|X(m)|^2) = 2*10*log10(|X(m)|)
%
% @param   X      The DFT frequency domain sequence.
% @param   norm   Optional parameter to normalize (true) the data or not
%                 (false, default).  Normalization is |X(m)|/|X(0)| so that
%                 the first value in the normalizesd log magnitude sequence
%                 is equal to 0dB.
%               
%
% @return  The power spectrum sequence in decibels.
%
% @throw   Error if the input sequence is not a vector.
%
% [1] Lyons, Richard G., Understanding Digital Signal Processing 2nd 
%     Edition, Prentice Hall 2004
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
function [P] = todB(X,norm)

    if ~exist('norm','var')
        norm = false;
    end
    
    Xmag = abs(X);
    
    if norm
        P = 20*log10(Xmag./max(Xmag));
    else
        P = 20*log10(Xmag);
    end

end