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
% http://www.movable-type.co.uk/scripts/latlong-vincenty.html
% http://www.sosmath.com/trig/Trig5/trig5/trig5.html
% http://williams.best.vwh.net/avform.htm
% http://www.ga.gov.au/earth-monitoring/geodesy.html
%
% [A] Direct and Inverse Solutions of Geodesics on the Ellipsoid with
%     Application of Nested Equations, T. Vincenty, April 1975 (Survey
%     Review XXII)
% [B] Geodetic Inverse Solution Between Antipodal Points, T. Vincenty,
%     August 1975.
%
% INPUTS:
% lat1     - latitude of first position in radians (positive north)
% lon1     - longitude of first position in radians (positive east)
% lat2     - latitude of second position in radians (positive north)
% lon2     - longitude of second position in radians (positive east)
%
% OUTPUTS:
% range    - range of second position from first in metres
% azimuth1 - azimuth of second position from first in radians (positive
%            north); is the initial bearing or forward azimuth
% azimuth2 - final bearing in radians (in direction from the first point to
%            the second point)
% err      - error code, successful (1) or failure (0)
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
function [ range, azimuth1, azimuth2, err ] = getRangeAzimuth( lat1, lon1, lat2, lon2 )

% load constants
wgs84_constants();

global SEMI_MAJOR_AXIS;
global SEMI_MINOR_AXIS;

a = SEMI_MAJOR_AXIS;
b = SEMI_MINOR_AXIS;

MAX_ITERATIONS = 100;

f = (a-b)/a; % flattening, inverse of INVERSE_FLATTENING

phi1   = lat1;
phi2   = lat2;

% initialize the outputs
range    = 0.0;
azimuth1 = 0.0;
azimuth2 = 0.0;
err      = 0; % 0 for failure

% capture the case where the two geogpraphic positions are coincident
if abs(lat1-lat2)<eps && abs(lon1-lon2)<eps
    % the positions are considered identical so the range and azimuth(s)
    % are zero
    err = 1;
end

% calculate the reduced latitude (U1)
tanU1 = (1.0 - f)*tan(phi1);
cosU1 = 1/sqrt(1 + tanU1^2); % due to 1 + tan^2x = sec^2x where secx = 1/cosx
sinU1 = tanU1*cosU1;
U1    = atan2(sinU1,cosU1);

% calculate the reduced latitude (U2)
tanU2 = (1.0 - f)*tan(phi2);
cosU2 = 1/sqrt(1 + tanU2^2); % due to 1 + tan^2x = sec^2x where secx = 1/cosx
sinU2 = tanU2*cosU2;
U2    = atan2(sinU2,cosU2);

% calculate the difference in longitude, positive east
L = lon2 - lon1;

if L > pi
    L = L - 2.0*pi;
elseif L < -pi
	L = L + 2.0*pi;
end

    
lambda  = L;                                                                % [13]
lambdaP = 2*pi; 

% iterate until little change in lambda
epsilon = 10.0^-12; 

% force at least one iteration of the loop
iter = 0;

    if abs(lambda) > pi

        % anti-podal points must be handled separately
        if L > 0
            LP = 2.0*pi - L;
        elseif L < 0
            LP = -2.0*pi - L;
        end

        lambdaP    = 0;
        cosSqAlpha = 0.5;
        cos2sigmaM = 0.0;
        sigma      = pi - abs(U1+U2);
        sinAlpha   = sqrt(1-cosSqAlpha);
        sinAlphaP  = pi;

        epsilon = 10.0^-12;

        % force at least one iteration of the loop
        iter2 = 0;

        while abs(sinAlpha-sinAlphaP)>epsilon || iter2 == 0

            iter2 = iter2 + 1;
            if iter2 > MAX_ITERATIONS
                fprintf('!!! Max iterations reached - no convergence !!!\n');
                break;
            end

            C = (1/16.0)*f*cosSqAlpha*(4.0+f*(4.0-3.0*cosSqAlpha));

            if cosSqAlpha == 0.0
                return;
            end

            cos2sigmaM = cos(sigma) - 2.0*sinU1*sinU2/cosSqAlpha;

            D = (1.0-C)*f*(sigma+C*sin(sigma)*(cos2sigmaM + C*cos(sigma)*   ...
                (-1.0+2.0*cos2sigmaM^2)));

            if D==0.0
                return;
            end

            sinAlphaP  = sinAlpha;
            sinAlpha   = (LP-lambdaP)/D;
            
            if abs(sinAlpha) > 1.0
                if sinAlpha < 0
                    sinAlpha = -1.0;
                else
                    sinAlpha = 1.0;
                end
            end
            
            cosSqAlpha = 1.0 - sinAlpha^2;

            denom = cosU1*cosU2;

            if denom == 0.0
                return;
            end

            sinLambdaP = sinAlpha*sin(sigma)/denom;
            cosLambdaP = sqrt(1.0 - sinLambdaP^2);

            lambdaP    = atan2(sinLambdaP,cosLambdaP);

            sinSqSigma = (cosU2*sin(lambdaP))^2 +                           ...
                (cosU1*sinU2+sinU1*cosU2*cos(lambdaP))^2;

            cosSqSigma = 1.0 - sinSqSigma;

            sigma = atan2(sqrt(sinSqSigma), sqrt(cosSqSigma));

        end

        sinAlpha1 = sinAlpha/cosU1;
        cosAlpha1 = sqrt(1.0-sinAlpha1^2);
        
        if (cosU1*sinU2+sinU1*cosU2*cosLambdaP) < 0
            cosAlpha1 = -cosAlpha1;
        end

        azimuth1 = atan2(sinAlpha1,cosAlpha1);

        denom = -sinU1*sin(sigma) + cosU1*cos(sigma)*cos(azimuth1);

        if abs(denom) == 0.0
            return;
        end

        azimuth2 = atan2(sinAlpha,denom);   
        
        % now calculate the distance
        e = (a^2 - b^2)/b^2;
        E = sqrt(1.0 + e*cosSqAlpha);
        F = (E-1.0)/(E+1.0);
        A = (1.0+0.25*F^2)/(1.0-F);
        B = F*(1.0-3/8*F^2);
        delSigma = B*sin(sigma)*(cos2sigmaM+0.25*B*(cos(sigma)*(-1.0+2.0*cos2sigmaM^2) - ...
            1/6*B*cos2sigmaM*(-3.0+4.0*sin(sigma)^2)*(-3.0+4.0*cos2sigmaM^2)));
        s = (1.0-f)*a*A*(sigma-delSigma);

    else

        while abs(lambda - lambdaP) > epsilon || iter == 0

            iter = iter + 1;
            if iter > MAX_ITERATIONS
                %fprintf('!!! Max iterations reached - no convergence !!!\n');
                break;
            end        


            sinLambda = sin(lambda);
            cosLambda = cos(lambda);

            sinSqSigma = (cosU2*sinLambda)^2 + (cosU1*sinU2 -                 ...
                sinU1*cosU2*cosLambda)^2;                                           % [14]

            sinSigma = sqrt(sinSqSigma);

            cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda;                         % [15]

            if abs(cosSigma) == 0.0
                return;
            end

            sigma = atan2(sinSigma, cosSigma);

            tanSigma = sinSigma/cosSigma;                                           %#ok<NASGU> % [16]

            if abs(sinSigma) == 0.0
                return;
            end

            sinAlpha = cosU1*cosU2*sinLambda/sinSigma;                              % [17]

            cosAlpha = sqrt(1.0 - sinAlpha^2);

            if abs(cosAlpha) == 0.0
                return;
            end

            cos2sigmaM = cosSigma - 2.0*sinU1*sinU2/cosAlpha^2;                     % [18]

            % store the previous result
            lambdaP = lambda;

            C = (f/16.0)*cosAlpha^2*(4.0+f*(4.0-3.0*cosAlpha^2));                   % [10]
            lambda = L + (1.0-C)*f*sinAlpha*(sigma+C*sinSigma*(cos2sigmaM+      ... % [11]
                cosSigma*C*(-1.0+2.0*cos2sigmaM^2)));  

            sigma = atan2(sinSigma, cosSigma);   
        end

        uSqrd = cosAlpha^2 * (a^2 - b^2)/b^2;

        A = 1.0+(uSqrd/16384.0)*(4096.0+uSqrd*(-768.0+uSqrd*(320.0-175.0*uSqrd)));  % [3]
        B = (uSqrd/1024.0)*(256.0+uSqrd*(-128.0 + uSqrd*(74.0-47.0*uSqrd)));        % [4]

        deltaSigma = B*sinSigma*(cos2sigmaM + B/4*(cosSigma*                    ...
                (-1.0+2.0*cos2sigmaM^2) - B/6*cos2sigmaM*(-3.0+4.0*             ...
                sinSigma^2)*(-3.0+4.0*cos2sigmaM^2)));                              % [6]

        s = b*A*(sigma - deltaSigma);

        denom = cosU1*sinU2 - sinU1*cosU2*cos(lambda);

        azimuth1 = atan2(cosU2*sin(lambda),denom);                                  % [20]
        
        % for azimuth2 (from P2 back to P1) we set lambda to be the
        % negative of itself; with this in mind, sin(-lambda) =
        % -sin(lambda) and cos(-lambda) = cos(lambda); this results in
        % something a little different than [21] in the note
        denom = sinU1*cosU2 - cosU1*sinU2*cos(lambda);

        azimuth2 = atan2(-cosU1*sin(lambda),denom);                                  % [21]
    end
% successful range and azimuths
range    = s;
err      = 1;

end

