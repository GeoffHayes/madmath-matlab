% Function to write binary data files of observation and track history data
% to validate the implementation of the Linear Kalman Filter in java.
%
% Input:
%   saveFilesToDir   The directory to save the binary data files to.
%
% 2013-08-19    Geoff Hayes     Initial Release.
%
function javaTrackerDriver(saveFilesToDir)

% close all open files
    fclose all;
    
    % run the tracker on its default target parameters
    [xyObservations, rangeObservations, brgObservations, truth, ...
        xyTracks, rngTracks, brgTracks] = targetTracker();

    dataTypeDir = 'double';
    
    saveFilesToDir = [saveFilesToDir '/' dataTypeDir '/'];
    
    if ~isdir(saveFilesToDir)
        fprintf('The directory %s is invalid!!\n', saveFilesToDir);
        return;
    end   
    
    % write data to the LKF validation file
    fod = fopen([saveFilesToDir 'lkfValidation.bin'], 'wb');
    
    if fod
        
        % assume one observation for every track update
        for i=1:size(xyTracks,2)
            
            obs = xyObservations{i};
            trk = xyTracks{i};
            
            % write the observation data to file
            fwrite(fod, obs.initTimeSecs, 'double', 'b');
            [rows,cols] = size(obs.z);
            fwrite(fod,rows,'int16','b');
            fwrite(fod,cols,'int16','b');
            fwrite(fod, obs.z', 'double', 'b');
            
            [rows,cols] = size(obs.R);
            fwrite(fod,rows,'int16','b');
            fwrite(fod,cols,'int16','b');            
            fwrite(fod, obs.R', 'double', 'b');
            
            % write the track data to file
            fwrite(fod, trk.lastUpdateTimeSecs, 'double', 'b');
            
            [rows,cols] = size(trk.F);
            fwrite(fod,rows,'int16','b');
            fwrite(fod,cols,'int16','b');
            fwrite(fod,trk.F','double','b');
            
            [rows,cols] = size(trk.Q);
            fwrite(fod,rows,'int16','b');
            fwrite(fod,cols,'int16','b');
            fwrite(fod,trk.Q','double','b');
            
            [rows,cols] = size(trk.xp);
            fwrite(fod,rows,'int16','b');
            fwrite(fod,cols,'int16','b');
            
            if rows>0 && cols > 0
                fwrite(fod,trk.xp','double','b');
            end
            
            [rows,cols] = size(trk.Pp);
            fwrite(fod,rows,'int16','b');
            fwrite(fod,cols,'int16','b');
            if rows>0 && cols > 0
                fwrite(fod,trk.Pp','double','b');
            end
            
            [rows,cols] = size(trk.H);
            fwrite(fod,rows,'int16','b');
            fwrite(fod,cols,'int16','b');
            if rows>0 && cols > 0
                fwrite(fod,trk.H','double','b');
            end   
            
            [rows,cols] = size(trk.y);
            fwrite(fod,rows,'int16','b');
            fwrite(fod,cols,'int16','b');
            if rows>0 && cols > 0
                fwrite(fod,trk.y','double','b');
            end   
            
            [rows,cols] = size(trk.S);
            fwrite(fod,rows,'int16','b');
            fwrite(fod,cols,'int16','b');
            if rows>0 && cols > 0
                fwrite(fod,trk.S','double','b');
            end    
            
            [rows,cols] = size(trk.K);
            fwrite(fod,rows,'int16','b');
            fwrite(fod,cols,'int16','b');
            if rows>0 && cols > 0
                fwrite(fod,trk.K','double','b');
            end    
            
            [rows,cols] = size(trk.x);
            fwrite(fod,rows,'int16','b');
            fwrite(fod,cols,'int16','b');
            if rows>0 && cols > 0
                fwrite(fod,trk.x','double','b');
            end         
            
            [rows,cols] = size(trk.P);
            fwrite(fod,rows,'int16','b');
            fwrite(fod,cols,'int16','b');
            if rows>0 && cols > 0
                fwrite(fod,trk.P','double','b');
            end           
       
        end       
        
        fclose(fod);
    end
    
    