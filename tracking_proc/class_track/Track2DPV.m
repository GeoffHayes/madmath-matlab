%**************************************************************************
% Class encapsulating a 2D (Cartesian (x,y)) Position and Velocity 
% track.  An instance of this class is typically instantiated from a 
% single observation, and the track is then periodically or 
% aperiodically updated with future observations.
%**************************************************************************
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
classdef Track2DPV < Track
    
    properties
        
        x  = zeros(Track2DPV.DIMS,1);
        P  = zeros(Track2DPV.DIMS,Track2DPV.DIMS);
        xp = zeros(Track2DPV.DIMS,1);
        Pp = zeros(Track2DPV.DIMS,Track2DPV.DIMS);
        H  = [];
        y  = [];
        K  = [];
        S  = [];
        F  = zeros(Track2DPV.DIMS,Track2DPV.DIMS);
        Q  = zeros(Track2DPV.DIMS,Track2DPV.DIMS);
    end
    
    properties (Constant, Access=public)
       DIMS  = 4;
       XPOS  = 1;
       YPOS  = 2;
       VX    = 3;
       VY    = 4;
    end  
    
    methods
        % class constructor
        function [trk] = Track2DPV()
            
            % call the superclass constructor
            trk = trk@Track();
           
        end 
        
        function [trk] = init(trk, obs)
            
            if isa(obs,'Observation')
                % initialize the track with the observation
                trk = obs.initTrack(trk);
            elseif isa(obs,'Track2DPV')
                trk.x = obs.x;
                trk.P(Track2DPV.XPOS,Track2DPV.XPOS) = 500^2;
                trk.P(Track2DPV.YPOS,Track2DPV.YPOS) = 500^2; 
            end
            
            % set the track timestamps
            trk.initTimeSecs       = obs.initTimeSecs;
            trk.lastUpdateTimeSecs = obs.initTimeSecs;
            
            % initialize the covariance velocities
            trk.P(Track2DPV.VX,Track2DPV.VX) = 1.5^2;
            trk.P(Track2DPV.VY,Track2DPV.VY) = 1.5^2;
            
        end
        
        function [trk] = predict(trk, predictTimeSecs)
            
            % calculate the elapsed time
            elapsedTimeSecs = predictTimeSecs - trk.lastUpdateTimeSecs;
            
            F = eye(size(trk.P));
            F(Track2DPV.XPOS,Track2DPV.VX) = elapsedTimeSecs;
            F(Track2DPV.YPOS,Track2DPV.VY) = elapsedTimeSecs;

            dt = elapsedTimeSecs;
            Q = [dt^3/3 0      dt^2/2 0;
                 0      dt^3/3 0      dt^2/2;
                 dt^2/2 0      dt     0;
                 0      dt^2/2 0      dt];
             
            q = 0.16;
            
            trk.x = F*trk.x;
            trk.P = F*trk.P*F' + q*Q;
            
            trk.lastUpdateTimeSecs = predictTimeSecs;
            
            % save the predicted x, P and F and Q
            trk.xp = trk.x;
            trk.Pp = trk.P;
            trk.F  = F;
            trk.Q  = Q;

        end
        
        function [trk] = update(trk,obs)
            
            % get the mapping matrix H
            H = obs.getH(trk);
            
            % compute the innovation vector y=z-Hx
            y = obs.getY(trk,H);
            
            % compute the innovation covariance matrix S=HPH' + R
            S = H*trk.P*H' + obs.R;
            
            % compute the Kalman gain K=PH'inv(S)
            K = trk.P*H'*inv(S);
            
            % update the state vector x=x+Ky
            trk.x = trk.x + K*y;
            
            % update the covariance matrix P=(I-KH)P(I-KH)'+ KRK'
            I = eye(size(trk.P));
            trk.P = (I-K*H)*trk.P*(I-K*H)' + K*obs.R*K';
            
            % save the H, y, S and K
            trk.H = H;
            trk.y = y;
            trk.S = S;
            trk.K = K;
            
        end
        
        % output methods
        function [str] = asString(trk)
            str = [trk.asString@Track() ' x:' mat2str(trk.x) ' P:' mat2str(trk.P)];
        end
        
        function display(trk)
            disp(trk.asString()); 
        end
    end
end

