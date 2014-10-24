%**************************************************************************
% Abstract class to encapsulate common attributes (data members) for 
% tracks that are generally produced (predicted and updated) via any 
% filter (typically Kalman).
%**************************************************************************
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
classdef (Abstract) Track

    properties
        % Latitude and longitude position of the track.
        position;
        % Time at which the track is initialized (seconds).
        initTimeSecs;
        % Last time at which the track was predicted and/or updated
        % (seconds).
        lastUpdateTimeSecs;
    end
    
    properties (Abstract)
        % State vector.
        x;
        % Covariance matrix for the state vector.
        P;
        % predicted state vector
        xp;
        % predicated covariance matrix
        Pp;
        % mapping from track to observation space matrix
        H;
        % innovation or residual (difference between the observed and
        % actual/predicted state)
        y;
        % Kalman gain matrix
        K;
        % innovation covaraince matrix
        S;
        % transition matrix
        F;
        % process noise matrix
        Q;
    end
    
    methods (Abstract)
        init(trk, obs);
        predict(trk, predictTimeSecs);
        update(trk,obs);
    end
    
    methods
        % class constructor
        function [trk] = Track()

            trk.position           = ...
                LatLonPosition(0.0,0.0,LatLonPosition.RADIANS);
            trk.initTimeSecs       = 0.0;
            trk.lastUpdateTimeSecs = 0.0;
           
        end 
        
        % set methods
        function [trk] = setPostion(trk,pos)
            trk.position = pos;
        end
        
        function [trk] = setInitTime(trk,initTimeSecs)
            if initTimeSecs < 0.0
                initTimeSecs = 0.0;
            end
            trk.initTimeSecs = initTimeSecs;
        end
        
        function [trk] = setLastUpdateTime(lastUpdateTimeSecs)
            if lastUpdateTimeSecs < 0.0
                lastUpdateTimeSecs = 0.0;
            end
            trk.lastUpdateTimeSecs = lastUpdateTimeSecs;
        end
        
        % get methods
        function [pos] = getPosition(trk)
            pos = trk.position;
        end
        
        function [initTimeSecs] = getInitTimeSecs(trk)
            initTimeSecs = trk.initTimeSecs;
        end
        
        function [lastUpdateTimeSecs] = getLastUpdateTimeSecs(trk)
            lastUpdateTimeSecs = trk.lastUpdateTimeSecs;
        end
        
        % output methods
        function [str] = asString(trk)
            str = [trk.position.asString() ' ' num2str(trk.initTimeSecs) ' ' ...
                num2str(trk.lastUpdateTimeSecs)];
        end
        
        function display(trk)
            disp(trk.asString()); 
        end
    end  
end
