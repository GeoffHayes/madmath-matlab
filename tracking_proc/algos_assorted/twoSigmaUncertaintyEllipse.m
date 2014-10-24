%**************************************************************************
% Determines the two-sigma uncertainty ellipse of the object (track,
% target, observation, etc.) given its 2D position covariance matrix.
%
% @param   P   The 2D position covariance matrix with the position being a
%              (x,y) Cartesian coordinate relative to some origin.
%
% @return  The half-length of the ellipse along the x-axis (a), the
%          half-length the ellipse along the y-axis (b), and the 
%          orientation of the ellipse (angle of the ellipse semi-major axis 
%          relative to the x-axis).  Note that the half-length refers to
%          the length from the centre of the ellipse to the intersection of
%          the ellipse with the x-axis or y-axis.  The ellipse semi-major
%          axis is the axis (x or y) with the maximum of a and b.
%
% @throw   Error if no eigenvectors can be found.
%
% [1] Hary, Kenneth,Elementary Linear Algebra, Integer Press, 1995.
% [2] http://www.earth-time.org/projects/upb/public_docs/ErrorEllipses.pdf
% [3] http://math.stackexchange.com/questions/23596/why-is-the-eigenvector-of-
%         a-covariance-matrix-equal-to-a-principal-component
% [4] http://www.casaxps.com/help_manual/mathematics/
%         Matrices_and_EigenvectorsRev6.pdf
% [5] https://www.ngs.noaa.gov/PUBS_LIB/
%     AlgorithmsForConfidenceCirclesAndEllipses_TR_NOS107_CGS3.pdf
% [6] http://www.gpsinformation.net/main/errors.htm
% [7] http://users.tkk.fi/mvermeer/uncertainty.pdf
%
%**************************************************************************
% Code is written by author based upon noted references.  If you have any
% questions, comments, or find bugs, please feel free to email me at 
% geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
function [a,b,phi] = twoSigmaUncertaintyEllipse(P)
 
    % from [3], the "the eigenvector with the largest eigenvalue is the 
    % direction along which the data set has the maximum variance"...so we
    % need the eigenvalues and eigenvectors of the position covariance data
    
    % from [1]:
    
    % a vector x is an eigenvector if x is nonzero and satisfies 
    % R*x=lambda*x
    % ==> (R-lambda*I)x = 0 has non-trivial solutions (lambda~-0)
    % ==> (R-lambda*I) is not invertible
    % ==> det(R-lambda*I) = 0
    % P(lambda) = det(R-lambda*I) is called the characteristic polynomial
    % det(R-lambda*I) = 0 is called a characteristic equation

    % calculate B = R - lambda*I as:

    % b11 = sxx - lambda
    % b12 = sxy
    % b21 = syx
    % b22 = syy - lambda
    
    % calculate the determinant of B, det(B) as

    % |B| = b11*b22 - b12*b21
    %     = (sxx - lambda)*(syy - lambda) - sxy*syx
    %     = lambda^2 - (sxx+syy)lambda + (sxx*syy - sxy*syx)

    % use the quadratic formula to find the roots, r as

    % lambda = (-b +-sqrt(b^2-4ac))/2a where

    % ax^2 + bx + lambda = 0
    
    % extract the position data from the uncertainty ellipse - this
    % corresponds to the (position) random variables X and Y
    sxx = P(1,1);
    sxy = P(1,2);
    syy = P(2,2);

    a      = 1;
    b      = -(sxx+syy);
    lambda = sxx*syy - sxy*sxy;

    temp = b^2 - 4.0*a*lambda;

    if temp >= 0.0

        % initialize the vector of eigenvalues
        V = [(-b + sqrt(temp))/(2*a); (-b - sqrt(temp))/(2*a)]; 
        
        % are both eigenvalues positive
        if V(1)<=0 || V(2)<=0
            error('uncertainty_ellipse: at least one eigenvalue is non-zero');
        end

        % note that the eigenvectors for the two eigenvalues are
        % perpendicular to one another (due to the nature of eigenvectors)
        % and correspond then to the semi-major and semi-minor axes of the
        % ellipse
 
        % solve for (P-lambda*I)x = 0 using the first equation, x=[a b]'
        % (sxx-lambda)a + sxyb = 0
        
        % solve for the unit eigenvector (i.e. magnitude is one)
        % [a b]' ==> [cos(phi) sin(phi)]', phi is angle relative to x-axis
        
        % (sxx-lambds)(cos(phi)) + sxy(sin(phi)) = 0
        % ==> -(sxx-lambda)(cos(phi)) = sxy(sin(phi))
        % ==> -(sxx-lambda)/sxy       = sin(phi)/cos(phi)
        % ==>  (lambda-sxx)/sxy       = tan(phi)
        % ==>  ((((sxx+syy) +- sqrt((sxx+syy)^2 - 4(sxx*syy - sxy*sxy)))/2) - sxx)/sxy = tan(phi)
        % ==>  ((syy-sxx) +- sqrt((sxx-syy)^2 + 4(sxy*sxy)))/(2sxy) = tan(phi)
        % 
        % using tan trig identity of tan(2phi) = 2tan(phi)/(1-tan(phi)^2)
        %
        % ==> tan(2phi) = 2sxy/(sxx-syy);
        
        phi = 0.5*atan2(2*sxy,sxx-syy);
        
        % the semi-major axis is the square root of the larger of the two
        % eigenvalues; the semi-minor axis is the square root of the
        % smaller of the two eigenvalues
        
        V = sqrt(V);
        
        mx = max(V);
        mn = min(V);
        
        % assume that hte semi-major axis is always paralle to the x-axis
        a = mx*2.0; % multiply by two since 2-sigma 
                    % uncertainty ellipse (see [6])

        % semi-minor axis is parallel to the y-axis
        b = mn*2.0; % multiply by two since 2-sigma 
                    % uncertainty ellipse (see [6])      
      
    end    
end

