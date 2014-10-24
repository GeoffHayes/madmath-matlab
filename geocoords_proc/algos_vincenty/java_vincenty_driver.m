% Function to read in data files from the data directory to validate the
% Vincenty algorithms.
%**************************************************************************
% Code is written by author based upon noted references (if given), using 
% The MathWorks MATLAB function signature (if applicable) for convenience 
% only.
%
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************

function java_vincenty_driver()

    format long g
    
    EPSILON = 1e-11;

    % get a list of the text files from the data directory
    dataDir  = '/Users/geoff/Documents/MATLAB/position/data';
    dataToDir ='/Users/geoff/Development/java/vincenty/data/double/'
    fileList = dir(dataDir);
    
    testSet  = 0;
    
   % fod = fopen([dataToDir 'vincentyDataTest.bin'], 'wb');
    
%     if ~fod
%         fprintf('Could not open file %s\n', dataToDir);
%         return;
%     end
       
    for i=1:size(fileList,1)
       if ~isempty(strfind(fileList(i).name,'.txt'))
          % read the data from file - assumes all strings
          [lat, lon, city, country] = textread([dataDir '/' fileList(i).name], '%s%s%s%s', 'delimiter', '\t'); %#ok<REMFF1>
          
          fprintf('Validating data for file %s......|', fileList(i).name);
          numTests = 0;
          numPass  = 0;
          
          errorList = {};
          numErrors = 0;
          
          % now compute the range and bearing/azimuth between each pair of
          % cities
          for j=1:size(lat,1)
              
              testSet  = testSet+1;
              testCase = 0;
             
              latRadsJ = getLatitude(char(lat(j)));
              lonRadsJ = getLongitude(char(lon(j)));

              for k=1:size(lat,1)
                  
                  % update the console progress indicator
                  if mod(numTests,2) == 0
                      fprintf('\b-');
                  else
                      fprintf('\b|');
                  end
                  
                  if numTests==787
                      g=2;
                  end
                 
                  numTests = numTests + 1;
                  testCase = testCase + 1;

                  latRadsK = getLatitude(char(lat(k)));
                  lonRadsK = getLongitude(char(lon(k)));
                  
                  % calculate the range and azimuth
                  [range, azimuth, azimuth2, err1] = getRangeAzimuth(...
                      latRadsJ, lonRadsJ, latRadsK, lonRadsK); %#ok<ASGLU>
                  
                  % write to file
                  data = [latRadsJ lonRadsJ latRadsK lonRadsK range azimuth azimuth2];
                 % fwrite(fod,data,'double','b');
                  
                  % calculate the lat and lon using the range and azimuth
                  [latRads, lonRads, newAzimuth, err2] = getLatLong(...
                      latRadsJ, lonRadsJ, range, azimuth); %#ok<ASGLU>

                  % check to see if all calculations were performed
                  % correctly
                  if ~err1 || ~err2 || abs(latRads-latRadsK) > EPSILON || ...
                          abs(lonRads-lonRadsK) > EPSILON
                      
                      numErrors = numErrors + 1;
                      errorList{numErrors}.testSet  = testSet; %#ok<*AGROW>
                      errorList{numErrors}.testCase = testCase;
                      errorList{numErrors}.latRef   = char(lat(j));
                      errorList{numErrors}.lonRef   = char(lon(j));
                      errorList{numErrors}.latDest  = char(lat(k));
                      errorList{numErrors}.lonDest  = char(lon(k));
                      errorList{numErrors}.cityRef  = char(city(j));
                      errorList{numErrors}.ctryRef  = char(country(j));
                      errorList{numErrors}.cityDest = char(city(k));
                      errorList{numErrors}.ctryDest = char(country(k));
                      
                      if abs(latRads-latRadsK) > EPSILON 
                          errorList{numErrors}.errorType = 'lat diff';
                      elseif abs(lonRads-lonRadsK) > EPSILON
                          errorList{numErrors}.errorType = 'lon diff';
                      elseif ~err1
                          errorList{numErrors}.errorType = 'getRangeAzimuth error';
                      elseif ~err2
                          errorList{numErrors}.errorType = 'getLatLon error';
                      end

                  else
                      numPass = numPass + 1;
                  end  
              end
          end
          
          passRate = numPass/numTests;
          
          fprintf('\b%d/%d  (%.4f%%)\n', numPass, numTests, passRate*100.0);
          
          % write out the errors
          fprintf('\nErrors for this file:\n');
          for m=1:size(errorList,2)
             fprintf('%d  TS=%d TC=%d Ref:%s %s (%s, %s) Dest: %s %s (%s, %s)  Error: %s\n', ...  
                      m, errorList{m}.testSet, errorList{m}.testCase, errorList{m}.latRef,   ...
                      errorList{m}.lonRef, errorList{m}.cityRef, errorList{m}.ctryRef,       ...
                      errorList{m}.latDest, errorList{m}.lonDest, errorList{m}.cityDest,     ...
                      errorList{m}.ctryDest, errorList{m}.errorType);
          end
       end
    end
    
    fclose(fod);
end

% Function to convert the latitude string into the equivalent in radians.
function [latRads] = getLatitude(latStr)

    DEGTORAD = pi/180.0;

    [degs, remain] = strtok(latStr, '.');
    [mins, remain] = strtok(remain, '.');
    [dir,  remain] = strtok(remain, '.'); %#ok<NASGU>
    
    % convert the degrees and minutes (and seconds if available) into the
    % radian equivalent
    latRads = str2double(degs)*DEGTORAD + str2double(mins)/60.0*DEGTORAD;
    
    % apply the direction
    if strcmpi(dir,'S')
        latRads = -latRads;
    end

end

% Function to convert the longitude string into the equivalent in radians.
function [lonRads] = getLongitude(lonStr)

    DEGTORAD = pi/180.0;

    [degs, remain] = strtok(lonStr, '.');
    [mins, remain] = strtok(remain, '.');
    [dir,  remain] = strtok(remain, '.'); %#ok<NASGU>
    
    % convert the degrees and minutes (and seconds if available) into the
    % radian equivalent
    lonRads = str2double(degs)*DEGTORAD + str2double(mins)/60.0*DEGTORAD;
    
    % apply the direction
    if strcmpi(dir,'W')
        lonRads = -lonRads;
    end

end
