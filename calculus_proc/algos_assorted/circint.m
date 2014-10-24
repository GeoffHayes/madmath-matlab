%CIRCINT Calculates the two intersection points (if they exist) for two
%circles.
%   CIRCINT(x1,y1,r1,x2,y2,r2), calculates the intersection of two circles
%   centred at (x1,y1) and (x2,y2) with radii of r1 and r2 respectively.
%   If no intersection exists, an error is thrown.
%
%   [ix,iy]=CIRCINT(x1,y1,r1,x2,y2,r2), calculates the intersection of two
%   circles centred at (x1,y1) and (x2,y2) with radii 
%   of r1 and r2 respectively.  The x- and y-intersection points are
%   returned (if they exist).
%
%   Examples:
%       [ix1,iy1,ix2,iy2]=CIRCINT(0,0,5,5,5,10);
%
%   If you have any questions, comments, or find bugs, please feel free to 
%   email me at geoff.hayes74@gmail.com.
%
%   Geoff Hayes 2014
function [ix,iy] = circint(x1,y1,r1,x2,y2,r2)

    % if same centre, then they can't intersect unless they have the
    % same radius
    if ~(x1==x2 && y1==y2)

        u = r1^2-r2^2 - (x1^2-x2^2) - (y1^2-y2^2);

        % expand the equations of the circle for each circles
        % subtract one from the other and solve for y or x depending
        % upon whether y1==y2 or x1==x2

        if y1~=y2
            % solve for x in either equation
            v = -(x1-x2)/(y1-y2);
            u = u/(-2*(y1-y2));
            a = 1+v^2;
            b = (2*u*v - 2*y1*v - 2*x1);
            c = x1^2 + u^2 - 2*y1*u + y1^2 - r1^2;
            ix1 = (-b + sqrt(b^2-4*a*c))/(2*a);
            ix2 = (-b - sqrt(b^2-4*a*c))/(2*a);

            if ix2<ix1
                t=ix1;
                ix1=ix2;
                ix2=t;
            end

            if isreal(ix1) && isreal(ix2)
                iy1 = u + (ix1)*v;
                iy2 = u + (ix2)*v;
            else
                error('crcint - no intersection');
            end
            
            ix = [ix1 ix2];
            iy = [iy1 iy2];

        else
            % solve for y in either equation
            v = -(y1-y2)/(x1-x2); 
            u = u/(-2*(x1-x2));
            a = 1+v^2;
            b = (2*u*v - 2*x1*v - 2*y1);
            c = y1^2 + u^2 - 2*x1*u + x1^2 - r1^2;
            iy1 = (-b + sqrt(b^2-4*a*c))/(2*a);
            iy2 = (-b - sqrt(b^2-4*a*c))/(2*a);

            if isreal(iy1) && isreal(iy2)

                ix1 = u + (iy1)*v;
                ix2 = u + (iy2)*v;

                if ix2<ix1
                    t=ix1;
                    ix1=ix2;
                    ix2=t;
                    t=iy1;
                    iy1=iy2;
                    iy2=t;
                end 
            else
                error('crcint - no intersection');
            end
            
            ix = [ix1 ix2];
            iy = [iy1 iy2];
            
        end
    end   
end