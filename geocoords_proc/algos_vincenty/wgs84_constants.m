%**************************************************************************
% The new World Geodetic System is called WGS 84. It is currently the 
% reference system being used by the Global Positioning System. It is 
% geocentric and globally consistent within ±1 m. Current geodetic 
% realizations of the geocentric reference system family International 
% Terrestrial Reference System (ITRS) maintained by the IERS are 
% geocentric, and internally consistent, at the few-cm level, while still 
% being metre-level consistent with WGS 84.
%
% http://en.wikipedia.org/wiki/WGS-84
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
function wgs84_constants()

% The WGS 84 datum surface is an oblate spheroid (ellipsoid) with major 
% (transverse) radius a = 6378137 m at the equator and flattening f = 
% 1/298.257223563. The polar semi-minor (conjugate) radius b then 
% equals a times (1?f), or 6356752.3142 m.

% Flattening is a measure of the compression of a circle or sphere along a 
% diameter to form an ellipse or an ellipsoid of revolution (spheroid
% respectively. Other terms used are ellipticity, or oblateness. The usual 
% notation for flattening is f and its definition in terms of the semi-axes 
% of the resulting ellipse or ellipsoid is flattening = f = (a-b)/a.  The 
% compression factor is b/a in each case. For the ellipse, this factor is 
% also the aspect ratio of the ellipse.

global SEMI_MAJOR_AXIS;
global SEMI_MINOR_AXIS;
global INVERSE_FLATTENING;

SEMI_MAJOR_AXIS    = 6378137.0;      % metres
SEMI_MINOR_AXIS    = 6356752.314245; % metres
INVERSE_FLATTENING = 298.257223563;