% Function that simulates a targetTruth and creates xyObservations based on the
% movement of said targetTruth.  These xyObservations are then used to instantiate
% and update a track periodically.
function [xyObservations, rngObservations, brgObservations, targetTruth, xyTracks, rngTracks, brgTracks] = targetTracker()
    
    close all;
    
    randn('state', 112358);
    rand('state', 132135);
    
    % generate some targetTruth data
    [targetTruth, sensorTruth] = generateTruth();

    % observe the targetTruth
    [xyObservations, rngObservations, brgObservations, targetTruth] = observeTarget(targetTruth, sensorTruth);
   
    % track the targetTruth
    [xyTracks]  = trackTarget(xyObservations,  targetTruth, 1);
    [rngTracks] = trackTarget(rngObservations, targetTruth, 2);
    [brgTracks] = trackTarget(brgObservations, targetTruth, 3);
    
    [tracks] = trackTargetWithAllData(rngObservations, brgObservations, targetTruth, 4);

end

% Function to follow a targetTruth, generating targetTruth and xyObservations as the
% targetTruth transits.
function [xyObservations, rngObservations, brgObservations, targetTruth] = observeTarget(targetTruth, sensorTruth)
    
    % uncertainties
    posUncertaintyMtrs   = 500.0;
    brgUncertaintyRads   = 1.0*pi/180.0;

    % observation and targetTruth arrays
    xyObservations    = {};
    rngObservations   = {};
    brgObservations   = {};

    % follow a straight line targetTruth, generating an observation every
    % timeDelta seconds
    for i=1:size(targetTruth,2)

        obs       = Observation2DP();
        rngObs    = ObservationRange();
        brgObs    = ObservationBearing();
        posNoisey = LatLonPosition;
        
        rangeMtrs      = 0.0;
        rangeMtrsNoisy = 0.0;
        
        if i==1
           
            posNoisey             = targetTruth{i}.pos;
            
            obs.position          = posNoisey;
            obs.initTimeSecs      = targetTruth{i}.timeSecs;

        else
            
            % get the range and azimuth of the target relative to its
            % previous position
            [rangeMtrs, aziFwdRads, aziRevRads, err] =                         ...
                getRangeAzimuth(                                               ...
                    targetTruth{i-1}.pos.getLatitude(LatLonPosition.RADIANS),  ...
                    targetTruth{i-1}.pos.getLongitude(LatLonPosition.RADIANS), ...
                    targetTruth{i}.pos.getLatitude(LatLonPosition.RADIANS),    ...
                    targetTruth{i}.pos.getLongitude(LatLonPosition.RADIANS));
            
            % add some noise
            rangeMtrsNoisy = rangeMtrs;
            %rangeMtrsNoisy = abs(rangeMtrs + posUncertaintyMtrs*randn());
            
            aziFwdRadsNoisy = aziFwdRads;
            %aziFwdRadsNoisy = aziFwdRads + brgUncertaintyRads*randn();
            
            % calculate the new position of the targetTruth with noise
            [lat, lon, revAzi, err] = getLatLong(                          ...
                targetTruth{i-1}.pos.getLatitude(LatLonPosition.RADIANS),  ...
                targetTruth{i-1}.pos.getLongitude(LatLonPosition.RADIANS), ...
                rangeMtrsNoisy,                                            ...
                aziFwdRadsNoisy);
            
            posNoisy = LatLonPosition(lat, lon, LatLonPosition.RADIANS);
            
            % create the 2D observation
            obs.position     = posNoisy;
            obs.initTimeSecs = targetTruth{i}.timeSecs;
            
            % calculate the new position of the targetTruth without noise
            [lat, lon, revAzi, err] = getLatLong(                          ...
                targetTruth{i-1}.pos.getLatitude(LatLonPosition.RADIANS),  ...
                targetTruth{i-1}.pos.getLongitude(LatLonPosition.RADIANS), ...
                rangeMtrs,                                                 ...
                aziFwdRads); 
        end
        
        % determine the state vector (z) and covaraince matrix (R) for the
        % observation
        obs.R = [posUncertaintyMtrs*posUncertaintyMtrs 0.0;       ...
                 0.0 posUncertaintyMtrs*posUncertaintyMtrs];
             
        [range, fwdAzi, revAzi, err] = getRangeAzimuth(           ...
            sensorTruth{1}.getLatitude(LatLonPosition.RADIANS),   ...
            sensorTruth{1}.getLongitude(LatLonPosition.RADIANS),  ...
            obs.position.getLatitude(LatLonPosition.RADIANS),     ...
            obs.position.getLongitude(LatLonPosition.RADIANS));
        
        % use the range and azimuth to determine the (x,y) coordinate
        x = range*sin(fwdAzi);
        y = range*cos(fwdAzi);
        
        obs.z = [x; y];
        
        xyObservations{i}   = obs;
        
        % find the range and fwdAzi from the origin (sensor 1) and the true
        % target position
        [range, fwdAzi, revAzi, err] = getRangeAzimuth(           ...
            sensorTruth{1}.getLatitude(LatLonPosition.RADIANS),   ...
            sensorTruth{1}.getLongitude(LatLonPosition.RADIANS),  ...
            targetTruth{i}.pos.getLatitude(LatLonPosition.RADIANS),       ...
            targetTruth{i}.pos.getLongitude(LatLonPosition.RADIANS));
        
        x = range*sin(fwdAzi);
        y = range*cos(fwdAzi);
        
        % update the range and bearing observations using the true range
        % and bearing respectively from the sensor to the position of the 
        % target without noise
        
        rngObs.sensorPos    = [0 0]'; 
        rngObs.initTimeSecs = targetTruth{i}.timeSecs;
        
        brgObs.sensorPos    = [0 0]'; 
        brgObs.initTimeSecs = targetTruth{i}.timeSecs;
        
        % add some noise to the range and bearing
        rngUncertaintyMtrs = 150.0;    
        brgUncertaintyRads = 1.0*pi/180.0;
        
        range = range;
        %range = abs(range + rngUncertaintyMtrs*randn());
        rngObs.z = range;
        rngObs.R = rngUncertaintyMtrs*rngUncertaintyMtrs;
        
        fwdAzi = fwdAzi;
        %fwdAzi = fwdAzi + brgUncertaintyRads*randn();
        brgObs.z = fwdAzi;
        brgObs.R = brgUncertaintyRads*brgUncertaintyRads;        
        
        rngObservations{2*i-1} = rngObs;
        brgObservations{2*i-1} = brgObs;
        
        % now add an observation from the same time but from a second
        % sensor
        [range, fwdAzi, revAzi, err] = getRangeAzimuth(               ...
            sensorTruth{2}.getLatitude(LatLonPosition.RADIANS),       ...
            sensorTruth{2}.getLongitude(LatLonPosition.RADIANS),      ...
            targetTruth{i}.pos.getLatitude(LatLonPosition.RADIANS),   ...
            targetTruth{i}.pos.getLongitude(LatLonPosition.RADIANS));  
        
        range = range;
        %range = abs(range + rngUncertaintyMtrs*randn());
        rngObs.z = range;
        rngObs.R = rngUncertaintyMtrs*rngUncertaintyMtrs;
        
        fwdAzi = fwdAzi;
        %fwdAzi = fwdAzi + brgUncertaintyRads*randn();
        brgObs.z = fwdAzi;
        brgObs.R = brgUncertaintyRads*brgUncertaintyRads;     
        
        % need to calculate the position of the second sensor relative to
        % the first
        [range, fwdAzi, revAzi, err] = getRangeAzimuth(            ...
            sensorTruth{1}.getLatitude(LatLonPosition.RADIANS),    ...
            sensorTruth{1}.getLongitude(LatLonPosition.RADIANS),   ...
            sensorTruth{2}.getLatitude(LatLonPosition.RADIANS),    ...
            sensorTruth{2}.getLongitude(LatLonPosition.RADIANS));  
        
        xs = range*sin(fwdAzi);
        ys = range*cos(fwdAzi);
        
        rngObs.sensorPos    = [xs ys]'; 
        brgObs.sensorPos    = [xs ys]';
        
        rngObservations{2*i} = rngObs;
        brgObservations{2*i} = brgObs;
        
    end

end

function [tracks] = trackTarget(observations, targetTruth, id)

    tracks = {};
    
    trk = Track2DPV();
    
    for i=1:size(observations,2)

        obs = observations{i};
        
        if i==1
            % initialize the track with the observation depending upon the
            % observation type
            if isa(obs,'Observation2DP')
                trk = trk.init(obs);
            else
                % cheat and initialize the track with the initial location
                % and speeds of the target
                temp              = Track2DPV();
                temp.x            = targetTruth{i}.state;
                temp.initTimeSecs = obs.initTimeSecs;
                trk               = trk.init(temp);
                
                 % predict the track to the time of the observation
                trk = trk.predict(obs.initTimeSecs);
            
                % update the track with the observation
                trk = trk.update(obs);
            end
        else
            % predict the track to the time of the observation
            trk = trk.predict(obs.initTimeSecs);
            
            % update the track with the observation
            trk = trk.update(obs);
        end
        
        tracks{i} = trk;
        
        % output the speed and course
        speed  = sqrt(trk.x(3)^2 + trk.x(4)^2);
        course = atan2(trk.x(3), trk.x(4))*180.0/pi;
        
        fprintf('speed=%.6f  course=%.6f\n', speed, course);
        
    end
    
    
    % plot the data
    figure(id);
    hold on;
    set(gca,'Color','k');
    for i=1:size(observations,2)
       if i<=size(targetTruth,2)
            plot(targetTruth{i}.state(1), targetTruth{i}.state(2), 'xg');
       end
       if isa(observations{i},'Observation2DP')
            plot(observations{i}.z(1),observations{i}.z(2), 'xr');
       else
           plot(observations{i}.sensorPos(1), observations{i}.sensorPos(2), '^r');
       end
       plot(tracks{i}.x(1), tracks{i}.x(2), 'ob');
    end
    
end

% function [targetTruth, sensorTruth] = generateTruth()
% 
%     targetTruth = {};
%     sensorTruth = {};
% 
%     % sensor position
%     sensorTruth{1} = LatLonPosition(45.4133172, -75.6901769, LatLonPosition.DEGREES);
%     sensorTruth{2} = LatLonPosition(45.7157752352646, -74.6323790982129, LatLonPosition.DEGREES);
%     
%     % targetTruth origin/initial position
%     targOrigin = LatLonPosition(45.4111700, -75.6981200, LatLonPosition.DEGREES);
%     
%     % targetTruth parameters
%     targSpeedMtrsPerSec = 25.60;
%     targHeadingRads     = 90*pi/180.0;
%     
%     timeSecs            = 123456789.0;
%     timeDeltaSecs       = 30.0;
% 
%     % follow a transiting line target, generating truth data for the target
%     % and sensors
%     for i=1:80
% 
%         if i>30 && i<=40
%             targHeadingRads = targHeadingRads + 18.0*pi/180.0;
%         end
%        
%         
%         if i>70 && i<=80
%             targHeadingRads = targHeadingRads + 18.0*pi/180.0;
%         end
%         
%         if i>80
%             targSpeedMtrsPerSec = 45.0;
%         end
% 
% %         if mod(i,2)==0
% %             aziFwdRads = aziFwdRads + 6.0*pi/180.0;
% %         end
%         
%         if i==1
%            
%             targPos = targOrigin;
%             
%         else
%             
%             % calculate the elapsed time since previous targetTruth (or
%             % observation) record
%             elapsedTimeSecs = timeSecs - targetTruth{i-1}.timeSecs;
%                    
%             % calculate the range from the previous position
%             rangeMtrs        = elapsedTimeSecs * targSpeedMtrsPerSec;
%             
%             % calculate the position of the target
%             [lat, lon, revAzi, err] = getLatLong(...
%                 targetTruth{i-1}.pos.getLatitude(LatLonPosition.RADIANS),  ...
%                 targetTruth{i-1}.pos.getLongitude(LatLonPosition.RADIANS), ...
%                 rangeMtrs,                                                 ...
%                 targHeadingRads);
%             
%             targPos = LatLonPosition(lat, lon, LatLonPosition.RADIANS);
%             
%         end
%         
%         
%         % set the target truth data
%         targetTruth{i}.pos      = targPos;
%         targetTruth{i}.timeSecs = timeSecs;
%         
%         % calculate the (x,y) relative to the first sensor which will act
%         % as the origin of the gaming area
%         [range, fwdAzi, revAzi, err] = getRangeAzimuth(          ...
%             sensorTruth{1}.getLatitude(LatLonPosition.RADIANS),  ...
%             sensorTruth{1}.getLongitude(LatLonPosition.RADIANS), ...
%             targPos.getLatitude(LatLonPosition.RADIANS),         ...
%             targPos.getLongitude(LatLonPosition.RADIANS));
%         
%         x = range*sin(fwdAzi);
%         y = range*cos(fwdAzi);
%         
%         vx = targSpeedMtrsPerSec*sin(targHeadingRads);
%         vy = targSpeedMtrsPerSec*cos(targHeadingRads);
%         
%         targetTruth{i}.state = [x y vx vy]';  
%         
%         % increment the time
%         timeSecs = timeSecs + timeDeltaSecs;
%         
%     end
% 
% end


function [tracks] = trackTargetWithAllData(bearings, ranges, targetTruth, id)

    tracks = {};
    
    trk = Track2DPV();
    
    for i=1:size(bearings,2)
        
        for j=1:2

            if j==1
                obs = bearings{i};
            else
                obs = ranges{i};
            end

            if i==1 && isempty(tracks)
                % initialize the track with the observation depending upon the
                % observation type
                if isa(obs,'Observation2DP')
                    trk = trk.init(obs);
                else
                    % cheat and initialize the track with the initial location
                    % and speeds of the target
                    temp              = Track2DPV();
                    temp.x            = targetTruth{i}.state;
                    temp.initTimeSecs = obs.initTimeSecs;
                    trk               = trk.init(temp);

                     % predict the track to the time of the observation
                    trk = trk.predict(obs.initTimeSecs);

                    % update the track with the observation
                    trk = trk.update(obs);
                end
            else
                % predict the track to the time of the observation
                trk = trk.predict(obs.initTimeSecs);

                % update the track with the observation
                trk = trk.update(obs);
            end
            
        end
        
        tracks{i} = trk;
        
        % output the speed and course
        speed  = sqrt(trk.x(3)^2 + trk.x(4)^2);
        course = atan2(trk.x(3), trk.x(4))*180.0/pi;
        
        fprintf('speed=%.6f  course=%.6f\n', speed, course);
        
    end
    
    
    % plot the data
    figure(id);
    hold on;
    set(gca,'Color','k');
    for i=1:size(bearings,2)
       if i<=size(targetTruth,2)
            plot(targetTruth{i}.state(1), targetTruth{i}.state(2), 'xg');
       end
       if isa(bearings{i},'Observation2DP')
            plot(bearings{i}.z(1),bearings{i}.z(2), 'xr');
       else
           plot(bearings{i}.sensorPos(1), bearings{i}.sensorPos(2), '^r');
       end
       plot(tracks{i}.x(1), tracks{i}.x(2), 'ob');
    end
    
end

