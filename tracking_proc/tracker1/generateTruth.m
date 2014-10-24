function [targetTruth, sensorTruth] = generateTruth()

    saveDataTo = '/Users/geoff/Development/android/workspace/SelfTrackerApp/res/raw/selftracker_demoa.bin';
    targetTruth = {};
    sensorTruth = {};

    % sensor position
    sensorTruth{1} = LatLonPosition(45.4133172, -75.6901769, LatLonPosition.DEGREES);
    sensorTruth{2} = LatLonPosition(45.7157752352646, -74.6323790982129, LatLonPosition.DEGREES);
    
    % targetTruth origin/initial position
    targOrigin = LatLonPosition(45.4111700, -75.6981200, LatLonPosition.DEGREES);
    
    % targetTruth parameters
    targSpeedMtrsPerSec = 25.00;
    targHeadingRads     = 45*pi/180.0;
    
    timeSecs            = 123456789.0;
    timeDeltaSecs       = 20.0;

    % follow a transiting line target, generating truth data for the target
    % and sensors
    for i=1:160

        if i>30 && i<=50
            targHeadingRads = targHeadingRads + 13.5*pi/180.0;
        end
        
        if i==51
            targSpeedMtrsPerSec = targSpeedMtrsPerSec*1.25;
        end
       
        
        if i==58
            targSpeedMtrsPerSec = targSpeedMtrsPerSec/1.25;
        end
            
        
        if i>57 && i<=77
            targHeadingRads = targHeadingRads - 13.5*pi/180.0;
        end
        
        if i==86
            targSpeedMtrsPerSec = targSpeedMtrsPerSec/2;
        end

%         if mod(i,2)==0
%             aziFwdRads = aziFwdRads + 6.0*pi/180.0;
%         end
        
        if i==1
           
            targPos = targOrigin;
            
        else
            
            % calculate the elapsed time since previous targetTruth (or
            % observation) record
            elapsedTimeSecs = timeSecs - targetTruth{i-1}.timeSecs;
                   
            % calculate the range from the previous position
            rangeMtrs        = elapsedTimeSecs * targSpeedMtrsPerSec;
            
            % calculate the position of the target
            [lat, lon, revAzi, err] = getLatLong(...
                targetTruth{i-1}.pos.getLatitude(LatLonPosition.RADIANS),  ...
                targetTruth{i-1}.pos.getLongitude(LatLonPosition.RADIANS), ...
                rangeMtrs,                                                 ...
                targHeadingRads);
            
            targPos = LatLonPosition(lat, lon, LatLonPosition.RADIANS);
            
        end
        
        
        % set the target truth data
        targetTruth{i}.pos      = targPos;
        targetTruth{i}.timeSecs = timeSecs;
        
        % calculate the (x,y) relative to the first sensor which will act
        % as the origin of the gaming area
        [range, fwdAzi, revAzi, err] = getRangeAzimuth(          ...
            sensorTruth{1}.getLatitude(LatLonPosition.RADIANS),  ...
            sensorTruth{1}.getLongitude(LatLonPosition.RADIANS), ...
            targPos.getLatitude(LatLonPosition.RADIANS),         ...
            targPos.getLongitude(LatLonPosition.RADIANS));
        
        x = range*sin(fwdAzi);
        y = range*cos(fwdAzi);
        
        vx = targSpeedMtrsPerSec*sin(targHeadingRads);
        vy = targSpeedMtrsPerSec*cos(targHeadingRads);
        
        targetTruth{i}.state = [x y vx vy]';  
        
        % increment the time
        timeSecs = timeSecs + timeDeltaSecs;
        
    end
    
%     fod = fopen(saveDataTo, 'wb');
    
    figure (1);
    hold on;
    set(gca,'Color','k');
    for j=1:size(targetTruth,2)
        state = targetTruth{j}.state;
        plot(state(1),state(2),'xb');
        
        fprintf('lat=%.6f lon=%.6f\n',targetTruth{j}.pos.getLatitude(LatLonPosition.DEGREES),...
           targetTruth{j}.pos.getLongitude(LatLonPosition.DEGREES));
        
%         if fod > 0
%             
%             fwrite(fod,targetTruth{j}.pos.getLatitude(LatLonPosition.DEGREES),'double','b');
%             fwrite(fod,targetTruth{j}.pos.getLongitude(LatLonPosition.DEGREES),'double','b');
%             
%         end
        
    end
    
    fclose all;

end