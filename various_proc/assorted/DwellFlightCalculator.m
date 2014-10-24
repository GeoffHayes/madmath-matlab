% Calculates the dwell and flight times for key pressing.
function DwellFlightCalculator

    % password length
    pwdLength   = 3;

    % password attempts
    pwdAttempts = 3;

    % the dwell time is the absolute time difference between when the 
    % key is pressed and that same key is released
    dwellTimes  = zeros(pwdAttempts,pwdLength);
    
    % the flight time is the time from releasing a key and pressing a
    % subsequent key
    flightTimes = zeros(pwdAttempts,pwdLength);

    meanFlightTimes = zeros(1,pwdLength);
    meanDwellTimes  = zeros(1,pwdLength);

    h = figure; 
    set(h,'KeyPressFcn',@KeyDownCb,'KeyReleaseFcn',@KeyUpCb) ;

    tic;

    numKeyDownEvents = 0;
    numKeyUpEvents   = 0;
    atPwdAttemptNum  = 1;

    function KeyDownCb(~,evnt)
        
        if atPwdAttemptNum<=pwdAttempts
        
            numKeyDownEvents = numKeyDownEvents + 1;

            fprintf('(%d,%d) Key: %s\n',atPwdAttemptNum,numKeyDownEvents,evnt.Key);

            % if this is the first key down event, then there is no was no key
            % pressed before, so the flight time is zero
            if atPwdAttemptNum==1 && numKeyDownEvents==1
                flightTime = 0;
            else
                % else a key has been released so we can compute the flight
                % time
                flightTime = toc;
            end

            % record the flight time
            flightTimes(atPwdAttemptNum,numKeyDownEvents) = flightTime;

            % start the "timer" for the dwell time
            tic;
        end
    end

    function KeyUpCb(~,evnt)
        
        if atPwdAttemptNum<=pwdAttempts
        
            numKeyUpEvents = numKeyUpEvents + 1;

            fprintf('(%d,%d) Key: %s\n',atPwdAttemptNum,numKeyUpEvents,evnt.Key);

            % since a key up follows a key down then we can compute the dwell
            % time
            dwellTime = toc;

            % record the dwell time
            dwellTimes(atPwdAttemptNum,numKeyUpEvents) = dwellTime;

            % move to the next password attempt if the number of key ups is
            % the same size as the password length
            if numKeyUpEvents==pwdLength
                
                atPwdAttemptNum  = atPwdAttemptNum + 1;
                numKeyUpEvents   = 0;
                numKeyDownEvents = 0;

                % have we reached the maximum number of password attempts?
                if atPwdAttemptNum>pwdAttempts
                    % calculate the means
                    meanFlightTimes = mean(flightTimes);
                    meanDwellTimes  = mean(dwellTimes);
                    
                    % write the results to a file (overwriting an existing
                    % file if it exists)
                    fid = fopen('dwellFlightStats.txt','wt');
                    
                    if fid>0
                        fmtStr = [repmat('\t%.10f',1,pwdLength) '\n'];
                        fprintf(fid,'Dwell Times:\n');
                        fprintf(fid,fmtStr,dwellTimes');
                        fprintf(fid,'Flight Times:\n');
                        fprintf(fid,fmtStr,flightTimes');
                        fprintf(fid,'Mean Dwell Times:\n');
                        fprintf(fid,fmtStr,meanDwellTimes');
                        fprintf(fid,'Mean Flight Times:\n');
                        fprintf(fid,fmtStr,meanFlightTimes');
                        fclose(fid); 
                    end
                    
                    % display all results to the console
                    disp('Dwell times:');
                    disp(dwellTimes);

                    disp('Flight times:');
                    disp(flightTimes);

                    disp('Mean dwell times:');
                    disp(meanDwellTimes);

                    disp('Mean flight times:');
                    disp(meanFlightTimes);                    
                else
                    % else start the timer for the flight time
                    tic;                
                end
            end
        end
    end
end