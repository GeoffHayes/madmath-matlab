%**************************************************************************
% An instance of the LatLonPostion class represents a geographical
% position given by a latitude and longitude.
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
classdef LatLonPosition
    
    %
    % See http://www.purplemath.com/modules/radians.htm for conversion of
    % degree latitude or longitude to degrees, mintues, seconds.  There are
    % 60 minutes per degree, and 60 seconds per minute.
    %
    
    properties
        latitude  = 0.0;
        longitude = 0.0;
    end
    
    % not quite an enumeration, but similar in nature
    properties (Constant)
       RADIANS   = 0;
       DEGREES   = 1;
    end    
    
    properties (Constant, Access=private)
       RADTODEG  = 180.0/pi;
       DEGTORAD  = pi/180.0;
       MAXLATDEG = 90.0;
       MAXLONDEG = 180.0;
       MAXLATRAD = pi/2;
       MAXLONRAD = pi;
    end  
    
    methods       
        
        % class constructor
        function [obj] = LatLonPosition(lat, lon, type)
            
            if ~exist('lat', 'var')
                lat = 0.0;
            end
            
            if ~exist('lon', 'var')
                lon = 0.0;
            end
            
            if ~exist('type', 'var')
                type = LatLonPosition.DEGREES;
            end
            
            
            obj = obj.setLatitude(lat,type);
            obj = obj.setLongitude(lon,type);
           
        end
        
        % set methods
        function [obj] = setLatitude(obj, lat, type)
            
            if ~exist('lat', 'var')
                lat = 0.0;
            end
            
            if ~exist('type', 'var')
                type = LatLonPosition.DEGREES;
            end
            
            if type == LatLonPosition.DEGREES
                if lat > LatLonPosition.MAXLATDEG
                    lat = LatLonPosition.MAXLATDEG;
                elseif lat < -LatLonPosition.MAXLATDEG;
                    lat = -LatLonPosition.MAXLATDEG;
                end
                
                obj.latitude  = lat*LatLonPosition.DEGTORAD;  
                
            elseif type == LatLonPosition.RADIANS
                
                if lat > LatLonPosition.MAXLATRAD
                    lat = LatLonPosition.MAXLATRAD;
                elseif lat < -LatLonPosition.MAXLATRAD;
                    lat = -LatLonPosition.MAXLATRAD;
                end
                
                obj.latitude  = lat;
            
            end
            
        end
        
        function [obj] = setLongitude(obj, lon, type)
            
            if ~exist('lon', 'var')
                lon = 0.0;
            end
            
            if ~exist('type', 'var')
                type = LatLonPosition.DEGREES;
            end
            
            if type == LatLonPosition.DEGREES
                if lon > LatLonPosition.MAXLONDEG
                    lon = LatLonPosition.MAXLONDEG;
                elseif lon < -LatLonPosition.MAXLONDEG;
                    lon = -LatLonPosition.MAXLONDEG;
                end
                
                obj.longitude  = lon*LatLonPosition.DEGTORAD;  
                
            elseif type == LatLonPosition.RADIANS
                
                if lon > LatLonPosition.MAXLONRAD
                    lon = LatLonPosition.MAXLONRAD;
                elseif lon < -LatLonPosition.MAXLONRAD;
                    lon = -LatLonPosition.MAXLONRAD;
                end
                
                obj.longitude  = lon;
            
            end
            
        end
        
        % get methods
        function [lat] = getLatitude(obj, type)
            
            if ~exist('type', 'var')
                type = LatLonPosition.DEGREES;
            end
            
            if type == LatLonPosition.DEGREES
                
                lat  = obj.latitude*LatLonPosition.RADTODEG;  
                
            elseif type == LatLonPosition.RADIANS
                
                lat  = obj.latitude;
            
            end
            
        end
        
        function [lon] = getLongitude(obj, type)
            
            if ~exist('type', 'var')
                type = LatLonPosition.DEGREES;
            end
            
            if type == LatLonPosition.DEGREES
                
                lon  = obj.longitude*LatLonPosition.RADTODEG;  
                
            elseif type == LatLonPosition.RADIANS
                
                lon  = obj.longitude;
            
            end
            
        end
        
        function [str] = asString(obj)

            % convert the latitude and longitude into degress, minutes,
            % seconds
            
            % convert the latitude and longitude into degrees
            latDegs = obj.latitude*LatLonPosition.RADTODEG;
            lonDegs = obj.longitude*LatLonPosition.RADTODEG;
            
            % convert to positive
            latDirection = 'N';
            lonDirection = 'E';
            
            if latDegs ~= abs(latDegs)
                latDegs = abs(latDegs);
                latDirection = 'S';
            end
            
            if lonDegs ~= abs(lonDegs)
                lonDegs = abs(lonDegs);
                lonDirection = 'W';
            end
            
            latTemp = latDegs - fix(latDegs);
            latMins = latTemp*60.0;
            latTemp = latMins - fix(latMins);
            latSecs = latTemp*60.0;
            
            % convert to integers
            latSecs = round(latSecs);
            latMins = fix(latMins);
            latDegs = fix(latDegs);
            
            if latSecs == 60
                latSecs = 0;
                latMins = latMins + 1;
                
                if latMins == 60
                    latMins = 0;
                    latDegs = latDegs + 1;
                    
                    if latDegs > LatLonPosition.MAXLATDEG
                        latDegs = LatLonPosition.MAXLATDEG;
                    end
                end
            end                
            
            lonTemp = lonDegs - fix(lonDegs);
            lonMins = lonTemp*60.0;
            lonTemp = lonMins - fix(lonMins);
            lonSecs = lonTemp*60.0;
            
            % convert to integers
            lonSecs = round(lonSecs);
            lonMins = fix(lonMins);
            lonDegs = fix(lonDegs);
            
            if lonSecs == 60
                lonSecs = 0;
                lonMins = lonMins + 1;
                
                if lonMins == 60
                    lonMins = 0;
                    lonDegs = lonDegs + 1;
                    
                    if lonDegs > LatLonPosition.MAXLONDEG
                        lonDegs = LatLonPosition.MAXLONDEG;
                    end
                end
            end
            
            str = [num2str(floor(latDegs),'%02d') '.' num2str(floor(latMins),'%02d') '.' ...
                   num2str(floor(latSecs),'%02d') latDirection ' '                      ...
                   num2str(floor(lonDegs),'%03d') '.' num2str(floor(lonMins),'%02d') '.' ...
                   num2str(floor(lonSecs),'%02d') lonDirection]; 
        end
        
        function display(obj)
           
            disp(obj.asString());
            
        end

    end
    
end

