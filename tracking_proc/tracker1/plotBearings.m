function plotBearings(tracks)

    MAXNUMSECONDSDISPLAYED = 30.0;
    
    prevTrack = tracks{1};

    minXCoord = 35;
    maxXCoord = 758;
    minYCoord = 0;    
    maxYCoord = 80;
    
    axisMinBearingDegs = 0.0;
    axisMaxBearingDegs = 360.0;
    
    minTimestampSecs  = prevTrack.lastUpdateTimeSecs+30;
    maxTimestampSecs  = minTimestampSecs + MAXNUMSECONDSDISPLAYED;
    
    timeScaleFactor = (maxXCoord-minXCoord)/ ...
        (maxTimestampSecs - minTimestampSecs);
            
    brgScaleFactor = (maxYCoord-minYCoord)/ ...
        (axisMaxBearingDegs-axisMinBearingDegs);
    
    startX = minXCoord;
    startY = maxYCoord;
    
    canvasUpdated = false;
    
    for i=2:size(tracks,2)
       
        nextTrack = tracks{i};
        
        exitEarly = false;

        startTimestampSecs = prevTrack.lastUpdateTimeSecs;

        stopTimestampSecs  = nextTrack.lastUpdateTimeSecs;

        % ensure that the bearing is within the correct range
        bearingStartDegs = atan2(prevTrack.x(3),prevTrack.x(4))*180.0/pi;
        
        if bearingStartDegs < 0.0
            bearingStartDegs = bearingStartDegs + 360.0;
        end

        bearingStartDegs = max(bearingStartDegs,axisMinBearingDegs);
        bearingStartDegs = min(bearingStartDegs,axisMaxBearingDegs);

        if ~canvasUpdated
            startY = startY - brgScaleFactor*bearingStartDegs;
            startY = min(startY, maxYCoord);
            startY = max(startY, minYCoord);
        end

        bearingStopDegs = atan2(nextTrack.x(3),nextTrack.x(4))*180.0/pi;
        
        if bearingStopDegs < 0.0
            bearingStopDegs = bearingStartDegs + 360.0;
        end

        bearingStopDegs = max(bearingStopDegs,axisMinBearingDegs);
        bearingStopDegs = min(bearingStopDegs,axisMaxBearingDegs);

        elapsedTimeSecs = ...
                stopTimestampSecs - startTimestampSecs;
            
        stopX = 0.0;
        stopY = 0.0;

        % should we draw this?
        fitsWithinWindow = true;
        
        if (startTimestampSecs > maxTimestampSecs)
            % the starting time is greater than the max plot time
            % so exit early
            fitsWithinWindow = false;
            
            exitEarly        = true;
            
        elseif (startTimestampSecs >= minTimestampSecs && ...
                 stopTimestampSecs  >  maxTimestampSecs)
           
            slope = ...
                ((bearingStopDegs - bearingStartDegs)*brgScaleFactor)/...
                ((elapsedTimeSecs)*timeScaleFactor);
            
            % the start time fits within the plot, but the
            % end time is outside of the plot so we need to
            % adjust the elapsed time stamp accordingly
            adjElapsedTimeSecs = maxTimestampSecs - ...
                    startTimestampSecs;

            % the start coordinate is fine, but the stop
            % coordinate is now on the edge at maxTimestampSecs
            stopX = maxXCoord;

            stopY = startY - adjElapsedTimeSecs*timeScaleFactor*slope;

            stopY = max(stopY,  minYCoord);
            stopY = min(stopY,  maxYCoord);

            % can skip any subsequent bearings since we have
            % painted until the end of the plot
            exitEarly = true;

        elseif (startTimestampSecs >= minTimestampSecs && ...
                 stopTimestampSecs  <=  maxTimestampSecs)
             
            % both times fit within the plot, so no adjustments
            % to the elapsed time is needed
            stopX = startX + timeScaleFactor*elapsedTimeSecs;

            stopY = maxYCoord - bearingStopDegs*brgScaleFactor;

            stopY = min(stopY,  maxYCoord);
            stopY = max(stopY,  minYCoord);
            
        elseif (startTimestampSecs < minTimestampSecs)
            
            slope = ...
                ((bearingStopDegs - bearingStartDegs)*brgScaleFactor)/...
                ((elapsedTimeSecs)*timeScaleFactor);
            
            adjElapsedTimeSecs = minTimestampSecs - ...
                    startTimestampSecs;
            
            % the start coordinate is then at the minimum
            % timestamp (i.e. the left of the plot)
            startX = minXCoord;

            startY = maxYCoord - adjElapsedTimeSecs*timeScaleFactor*slope;

            startY = max(startY,  minYCoord);
            startY = min(startY,  maxYCoord);

            if (stopTimestampSecs >= minTimestampSecs && ...
                stopTimestampSecs >  maxTimestampSecs)
                % so the plot window is enclosed by the start
                % and end times, so the elapsed time is just
                % the length of the window
                adjElapsedTimeSecs = MAXNUMSECONDSDISPLAYED;

                stopX  = maxXCoord;

                stopY = startY - adjElapsedTimeSecs*timeScaleFactor*slope;

                stopY = max(stopY,  minYCoord);
                stopY = min(stopY,  maxYCoord);   

                % exit early since the complete window is painted
                exitEarly = true;
                
            elseif (stopTimestampSecs >= minTimestampSecs && ...
                     stopTimestampSecs <=  maxTimestampSecs)
                 
                % the end time fits within the window but
                % the elapsed time needs to be adjusted
                adjElapsedTimeSecs = stopTimestampSecs - ...
                        minTimestampSecs;

                stopX = startX + timeScaleFactor*adjElapsedTimeSecs;

                stopY = startY - adjElapsedTimeSecs*timeScaleFactor*slope;

                stopY = max(stopY,  minYCoord);
                stopY = min(stopY,  maxYCoord); 

            else

                % the start and end times are both outside of
                % the window, so cannot be used
                fitsWithinWindow = false;
            end

        end

        if (stopX >= maxXCoord)
        
            exitEarly = true;
            stopX     = maxXCoord;  
        end
        
        if (startX == stopX)
            fitsWithinWindow = false;
        end

        % draw the line only if it has been fitted for the
        % window
        if (fitsWithinWindow)

%             % ensure that the correct colour is being used
%             if (!inHistoryMode)
%             {
%                 dataPaint.setColor(Color.GREEN);
%             }
%             else
%             {
%                 dataPaint.setColor(Color.MAGENTA);
%             }
%             canvas.drawLine((float)startX, (float)startY, (float)stopX, (float)stopY, dataPaint);

            startX = stopX;
            startY = stopY;
            
            canvasUpdated=true;
        end        

        if (exitEarly)
        
            break;
        end



        % save the previous data for the next iteration
        prevTrack = nextTrack;
    end


end

