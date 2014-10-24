%**************************************************************************
% The class encapsulates information concerning a two dimensional position 
% observation.
%**************************************************************************
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
classdef Observation2DP < Observation

    properties (Constant, Access=public)
       % Number of dimensions in the 2D position observation vector.
       DIMS  = 2;
       % Index of the x position (of 2D position) in the observation
       % vector.
       XPOS  = 1;
       % Index of the y position (of 2D position) in the observation
       % vector.
       YPOS  = 2;
    end  
    
    properties
        z = zeros(Observation2DP.DIMS,1);
        R = zeros(Observation2DP.DIMS,Observation2DP.DIMS); 
    end
    
    methods

        % class constructor
        function [obs] = Observation2DP()
            
            % call the superclass constructor
            obs = obs@Observation();
           
        end         
        
        function [H] = getH(obs,trk)
            
            if isa(trk,'Track2DPV')
               
                H = zeros(Observation2DP.DIMS,Track2DPV.DIMS);
                H(Observation2DP.XPOS,Observation2DP.XPOS) = 1.0;
                H(Observation2DP.YPOS,Observation2DP.YPOS) = 1.0;
                
            end
            
        end
        
        function [y] = getY(obs,trk,H)
            
            if isa(trk,'Track2DPV')
                y = obs.z - H*trk.x;
            end
        end
        
        
        function [trk] = initTrack(obs,trk)

            if isa(trk,'Track2DPV')

                trk.x = [obs.z; zeros(Observation2DP.DIMS,1)];
                trk.P = [obs.R zeros(Observation2DP.DIMS, Observation2DP.DIMS); ...
                    zeros(Observation2DP.DIMS,Observation2DP.DIMS)              ...
                    zeros(Observation2DP.DIMS,Observation2DP.DIMS)];
            end  
        end
        
        % output methods
        function [str] = asString(obs)
            str = [obs.asString@Observation() ' z:' mat2str(obs.z) ' R:' mat2str(obs.R)];;
        end
        
        function display(obs)
            disp(obs.asString());
        end        
    end
end

