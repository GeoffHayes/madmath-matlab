%DRAWSPHERE Draws a sphere for a given radius and origin.
%   DRAWSPHERE(radius), draws a sphere with radius units around the origin
%   (0,0,0).  The default spacing of points in a horizontal slice (xy 
%   plane) is 5 degrees.  The spacing of points in a vertical slice (zy 
%   plane) is 5 degrees.
%
%   DRAWSPHERE(radius,xyDeg,zyDeg), draws a sphere with radius units around 
%   the origin (0,0,0).  The spacing of points in a horizontal slice (xy 
%   plane) is xyDeg degrees.  The spacing of points in a vertical slice (zy 
%   plane) is zyDeg degrees.
%
%   DRAWSPHERE(radius,xyDeg,zyDeg,[u v w]), draws a sphere with radius units 
%   around the origin (u,v,w).  The spacing of points in a horizontal slice 
%   (xy plane) is xyDeg degrees.  The spacing of points in a vertical slice 
%   (zy plane) is zyDeg degrees.
%
%   [X,Y,Z] = DRAWSPHERE(radius,...), draws a sphere with radius units.
%   Returns the X, Y, and Z coordinates on the surface of the sphere such
%   that (X-u).^2 + (Y-v).^2 + (Z-w).^2 = radius^2 for the origin (u,v,w).
%
%   Examples:
%       drawsphere(10)
%       drawsphere(10,10,2)
%       drawsphere(10,10,5,[2 3 4])
%       drawsphere(10,36,5)  % decagon in xy plane
%       drawsphere(10,45,5)  % octagon in xy plane
%       drawsphere(10,72,5)  % pentagon in xy plane
%       drawsphere(10,90,5)  % square in xy plane
%       drawsphere(10,120,5) % triangles in xy plane
%
%   If you have any questions, comments, or find bugs, please feel free to 
%   email me at geoff.hayes74@gmail.com.
%
%   Geoff Hayes 2014
function [X,Y,Z] = drawsphere(radius,xyDeg,zyDeg,origin)

    if ~exist('xyDeg','var')
        xyDeg = 5;
    else
        xyDeg = max(1,min(360,abs(xyDeg)));
    end
    
    if ~exist('zyDeg','var')
        zyDeg = 5;
    else
        zyDeg = max(1,min(90,abs(zyDeg)));
    end
    
    if ~exist('origin','var')
        origin = [0 0 0];
    else
        if size(origin,2)>=3
            origin = origin(:,1:3);
        else
            origin = [origin(:,1:end) zeros(0,1,3-size(origin,2))];
        end
    end

    radius = abs(radius(1));
    
    Rxy  = [cosd(xyDeg)  -sind(xyDeg); ...
            sind(xyDeg)   cosd(xyDeg)];
        
    Rz   = [cosd(zyDeg)  -sind(zyDeg); ...
            sind(zyDeg)   cosd(zyDeg)]; 
      
    numPatchesXY = ceil(360/xyDeg+1);
    numPatchesZY = ceil(90/zyDeg+1);
    
    X = zeros(numPatchesXY,numPatchesZY*2-1);
    Y = zeros(numPatchesXY,numPatchesZY*2-1);
    Z = zeros(numPatchesXY,numPatchesZY*2-1);
    
    V = zeros(2,1);
    
    z    = zeros(numPatchesZY,1);
    z(1) = radius;
    V(1) = radius;
    V(2) = 0;
    for m=2:numPatchesZY
        V = Rz*V;
        z(m) = V(1);
    end
    z(end) = 0;

    z = [-z(1:end-1); flipud(z)];

    for k=1:length(z)
        
        Z(:,k) = repmat(z(k),numPatchesXY,1);
        
        V(1) = 0;
        V(2) = sqrt(radius^2 - z(k)^2);
        
        X(1,k) = V(1);
        Y(1,k) = V(2);
        
        for m=2:numPatchesXY
            
            V = Rxy*V;
            X(m,k) = V(1);
            Y(m,k) = V(2);
        end
    end
    
    % adjust for the origin
    X = X + origin(1);
    Y = Y + origin(2);
    Z = Z + origin(3);

    surf(X,Y,Z);
    axis equal;