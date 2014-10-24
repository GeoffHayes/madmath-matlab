%**************************************************************************
% Abstract class to encapsulate common attributes (data members) for 
% observations that are used to update tracks via a (say Kalman) 
% filter.
%**************************************************************************
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
classdef Observation

    properties
        % Latitude and longitude position of the observation.
        position;        
        % Time at which the observation was initialized (seconds).
        initTimeSecs;
    end
    
    properties (Abstract)
        % Observation state vector.
        z;
        % Covariance matrix for the observation state vector.
        R;
    end
    
    methods (Abstract)
        % Returns the mapping matrix from the true (target/track) space to
        % the observed space.
        getH(obs,trk);
        % Returns the innovation (or residual) of the observed value (the z
        % in the observation) and the actual value (in the track)
        getY(obs,trk,H);
        % Initializes the track with the observation data.
        initTrack(obs,trk);
    end
    
    methods
       
        % class constructor
        function [obs] = Observation()
            obs.initTimeSecs       = 0.0;
            obs.position           = ...
                LatLonPosition(0.0,0.0,LatLonPosition.RADIANS);            
        end 
        
        % Initializes the observation with the observed state (z), its
        % covariance matrix (R) and the initialization time (seconds).
        function [obs] = init(obs,z,R,initTimeSecs)
            if initTimeSecs < 0.0
                initTimeSecs = 0.0;
            end
            obs.initTimeSecs = initTimeSecs;
            obs.z            = z;
            obs.R            = R;
        end
        
        % get methods
        function [initTimeSecs] = getInitTimeSecs(obs)
            initTimeSecs = obs.initTimeSecs;
        end
        
        % output methods
        function [str] = asString(obs)
            
            str = [obs.position.asString() ' ' num2str(obs.initTimeSecs)];

        end
        
        function display(obs)
            disp(obs.asString());
        end
    end
end

