%**************************************************************************
% Vincenty's formulae are two related iterative methods used in geodesy to 
% calculate the distance between two points on the surface of a spheroid, 
% developed by Thaddeus Vincenty (1975a) They are based on the assumption 
% that the figure of the Earth is an oblate spheroid, and hence are more 
% accurate than methods such as great-circle distance which assume a 
% spherical Earth.
%
% The first (direct) method computes the location of a point which is a 
% given distance and azimuth (direction) from another point. The second 
% (inverse) method computes the geographical distance and azimuth between 
% two given points. They have been widely used in geodesy because they are 
% accurate to within 0.5 mm (0.020?) on the Earth ellipsoid.
%
% http://en.wikipedia.org/wiki/Vincenty's_formulae
% http://www.movable-type.co.uk/scripts/latlong-vincenty-direct.html
% http://www.sosmath.com/trig/Trig5/trig5/trig5.html
% http://www.ga.gov.au/earth-monitoring/geodesy.html
%
% [A] Direct and Inverse Solutions of Geodesics on the Ellipsoid with
%     Application of Nested Equations, T. Vincenty, April 1975 (Survey
%     Review XXII)
%
% INPUTS:
% refLat  - latitude of reference position in radians (positive north)
% refLon  - longitude of reference position in radians (positive east)
% range   - range from reference position to source in metres
% azimuth - angle from reference position to source in radians (positive
%           with respect to north)
%
% OUTPUTS:
% lat     - latitude of source in radians (positive north)
% lon     - longitude of source in radians (positive east)
% revAzi  - reverse azimuth in radians (positive north) from destination to
%           reference/source position
% err     - error code, successful (1) or failure (0)
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
function [lat, lon, revAzi, err] = getLatLong(refLat, refLon, range, azimuth)

% load constants
wgs84_constants();

global SEMI_MAJOR_AXIS;
global SEMI_MINOR_AXIS;

a = SEMI_MAJOR_AXIS;
b = SEMI_MINOR_AXIS;

MAX_ITERATIONS = 100;

f = (a-b)/a; % flattening, inverse of INVERSE_FLATTENING

s      = range;
alpha1 = azimuth;
phi1   = refLat;

% initialize the outputs
lat        = 0.0;
lon        = 0.0;
revAzi     = 0.0;
err        = 0; % 0 for failure

% return an error if the range is negative
if range < 0
    return;
elseif range < eps
    % the range is zero, so the lat and lon are identical to the reference
    % position
    lat = refLat;
    lon = refLon;   
    err = 1;
    return;
end

% calculate the reduced latitude (U1)
tanU1 = (1.0 - f)*tan(phi1);
cosU1 = 1/sqrt(1 + tanU1^2); % due to 1 + tan^2x = sec^2x where secx = 1/cosx
sinU1 = tanU1*cosU1;

% calcuate the angular distance on sphere from equator to reference
% position (sig1 of [A])
cosAlpha1 = cos(alpha1);
sinAlpha1 = sin(alpha1);
%tanSigma1 = tanU1/cosAlpha1;                                               % [1]
sigma1    = atan2(tanU1, cosAlpha1);

% calculate the azimuth of the geodesic at the equator
sinAlpha = cosU1*sinAlpha1;                                                 % [2]
cosAlpha = sqrt(1.0 - sinAlpha^2);

uSqrd = cosAlpha^2 * (a^2 - b^2)/b^2;

A = 1.0+(uSqrd/16384.0)*(4096.0+uSqrd*(-768.0+uSqrd*(320.0-175.0*uSqrd)));  % [3]
B = (uSqrd/1024.0)*(256.0+uSqrd*(-128.0 + uSqrd*(74.0-47.0*uSqrd)));        % [4]

sigma  = s/(b*A);
sigmaP = 2.0*pi;  % previous value

% iterate until little change in sigma
epsilon = 10.0^-12; 

% force at least one iteration of the loop
iter = 0;

while abs(sigma - sigmaP) > epsilon  || iter == 0
    
    iter = iter + 1;
    if iter > MAX_ITERATIONS
        %fprintf('!!! Max iterations reached - no convergence !!!\n');
        break;
    end
        
    % calculate the angular distance on the sphere from the equator to the
    % midpoint of the line (sigmaM)
    cos2sigmaM = cos(2.0*sigma1 + sigma);                                   % [5]
    
    sinSigma = sin(sigma);
    cosSigma = cos(sigma);
    
    deltaSigma = B*sinSigma*(cos2sigmaM + B/4*(cosSigma*         ...
        (-1.0+2.0*cos2sigmaM^2) - B/6*cos2sigmaM*(-3.0+4.0*      ...
        sinSigma^2)*(-3.0+4.0*cos2sigmaM^2)));                              % [6]
    
    sigmaP = sigma;
    
    sigma = s/(b*A) + deltaSigma;
    
end

sinSigma = sin(sigma);
cosSigma = cos(sigma);

% calculate the numer(ator) and denom(inator) of the formula that is used
% to determine the latitude (phi2)
denom = sinAlpha^2 + (sinU1*sinSigma - cosU1*cosSigma*cosAlpha1)^2;

% no solution if the denominator is less than or equal to zero
if denom <= 0.0
    return;
end

denom = sqrt(denom);

numer = sinU1*cosSigma + cosU1*sinSigma*cosAlpha1;

phi2 = atan2(numer, (1.0-f)*denom);                                         % [8]

% calculate the numer(ator) and denom(inator) of the formula that is used
% to determine the difference in longitude on an auxiliary sphere (lambda)
denom = cosU1*cosSigma - sinU1*sinSigma*cosAlpha1;

% no solution if the denominator is equal to zero
if denom == 0.0
    return;
end

numer = sinSigma*sinAlpha1;

lambda = atan2(numer, denom);                                               % [9]

% calculate the difference in longitude (L)
C = (f/16.0)*cosAlpha^2*(4.0+f*(4.0-3.0*cosAlpha^2));                       % [10]
L = lambda-(1.0-C)*f*sinAlpha*(sigma+C*sinSigma*(cos2sigmaM+            ... % [11]
    cosSigma*C*(-1.0+2.0*cos2sigmaM^2)));                            

% calculate the reverse azimuth (note that this is slightly different from
% the note, and the signs in the numerator and denominator have been
% reversed like they have been in the getRangeAzimuth function call - so
% [12] is different)
alpha2 = alpha1;
denom  = sinU1*sinSigma - cosU1*cosSigma*cosAlpha1;
if denom == 0.0
    % do nothing, leave the new azimuth as the original
else
    alpha2 = atan2(-sinAlpha, denom);                                       % [12]
end




% successful!
err        = 1;
lat        = phi2;
lon        = refLon + L;

if lon < -pi
    lon = lon + 2.0*pi;
elseif lon > pi
    lon = lon - 2.0*pi;
end

revAzi = alpha2;

end