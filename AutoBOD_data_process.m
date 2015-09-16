% AutoBOD_data_process.m
%
% Created 19 Aug 2015 by JRC under MATLAB R2015a; version history on GitHub
% Owner/author: James Collins, MIT-WHOI Joint Program, Woods Hole
% Oceanographic Institution, james.r.collins@aya.yale.edu
%
% Purpose: Read in, parse, and process dissolved oxygen data from the
% AutoBOD, a carousel developed in the Van Mooy Lab at WHOI to incubate,
% rotate, and make semi-continuous measurements of dissolved oxygen in
% multiple BOD bottles containing water samples
%
% Dependencies/required files:
%
% 1. An AutoBOD log and metadata file (.xlsx format with timestamps as ISO
% 8601) containing bottle sequence information, locations of individual
% data files, and other required inputs
%
% 2. The function linfit.m, available from
% http://www.whoi.edu/12.747/mfiles.html, or from the GitHub repo
% https://github.com/jamesrco/dependencies-useful-scripts
%
% Caveats:
%
% 1. Assumes data files (as generated by "capture text to file" function of
% your favorite terminal program) have the following consistent fixed-field
% format:
%
%   00033 000351 24.10 999.99 4 0 00033 001 0 08/13/15 18:07:17
%
%   Fields are (in order):
%       
%   1.  Amplitude (signal amplitude in uV)
%   2.  Phase (signal phase shift)
%   3.  Temperature (deg C)
%   4.  Oxygen concentration (for AutoBOD deployments prior to 9/10/15,
%       these values are in % saturation with the decimal two places from
%       the right termination of the string)
%
%       Units and decimal position are specified by the PreSens OEM unit
%       settings oxyu and ores. The values of these two settings at the
%       time of a given deployment should be recorded by user in the
%       deployment metadata spreadsheet.
%       
%       oxyu                      ores (for an example measurement "2341")
%       0*      % air sat         0     2341
%       1       % O2              1     234.1
%       2       hPa               2*    23.41
%       3       Torr              3     2.341
%       4       mg/L (ppm)        4     0.2341
%       5       micromoles per L
%
%       * OEM unit defaults
%
%   5.  "Error" (a "0" indicates that the unit is in range; '4' or
%       'amplitude too low' is a common error when light source does not
%       see a dot)
%   6.  "Bottle lock" (a "1" indicates AutoBOD is locked on a bottle and
%       not transiting)
%   7.  Amplitude again (this is a diagnostic field that; it is a
%       recalculated subroutine value used to to detect the edge of the
%       dot & should match the first field exactly)
%   8.  Bottle counter (not fully developed, as of 8/24/15; really at
%       this point just counts up to 255 and then starts over)
%   9.  Home detect (a "1" in this field indicates the carousel is
%       passing through the home sensor under the platter; this is how
%       we can keep track of the bottles for now)
%   10. Date (from computer clock)
%   11. Time (from computer clock)
%
% 2. Assumes user started AutoBOD with carousel at proper "home" position
% and that the first slug of "good" data appearing in the data file is for
% the bottle given in the log as bottle #1 for that deployment; if this is
% not the case, user may want to truncate the data file manually
%
% 3. Execution of script does not depend on the exact incubation start/end
% times given by the user in the deployment log, since the computer clock
% of the machine on which the data was recorded might have been different
% than the local time at the site; instead, the start/end times given in
% the log are used to calculate the length of the period of data from which
% results are calculated (where the first entry in the data file is t = 0)
%
% 4. The script includes several code snippets (user must select),
% depending on the format of the source DO data. Several scenarios exist:
%    
%   1. Onboard calculations were made as micromoles/L (OEM unit setting
%      (oxyu0005), so DO values can be imported directly
%   2. Onboard calculations were made as % air sat (oxyu0000, the default),
%      so conversion using concurrent temp and pressure measurements is
%      necessary
%   3. User desires calculation of values in micromoles/L directly from the
%      raw phase and amplitude measurements in the data file
%
% Note that in any of the above scenarios, values will still have to be
% adjusted for salinity (accomplished by a later section of code)

%% Clean things up and prep workspace

close all;
clear all;
clf;

%% User: specify some necessary file locations and parameters

% Name and location of bottle log & metadata file
AutoBOD_LogFile = '/Users/jrcollins/Dropbox/Cruises & projects/PHORCYS & AutoBOD/Data/AutoBOD_master_metadata.xlsx';

% Specify the deployments in the log file for which we want to read &
% process data; should be one or more of the deployments in the
% Deployment_metadata field "Deployment_ID"; if left unspecified, script
% will run through all deployments in the log file

%Desired_Deployments = ['Iselin_PHORCYS_2015_2'];

Max_RECD = 15; % The maximum allowable read error code deviations
              % this is the maximum number of bad reads allowable in a
              % stretch of "good" data (i.e., when locked onto a spot in a
              % particular bottle), or the number of erroneously good reads
              % in a "bad" (i.e., between bottles) stretch of data; necessary
              % because the measurement hardware isn't perfect
              
Num_bots = 12; % Number of bottles in the carousel (12)

% Specify some coefficients (to be used later on) for correcting DO values
% for salinity, from:
% Garcia HE, Gordon LI (1992), Oxygen solubility in seawater: Better
% fitting equations, Limnol. Oceanogr. 37:1307-1312

B_0 = -6.24097E-03;
B_1 = -6.93498E-03;
B_2 = -6.90358E-03;
B_3 = -4.29155E-03;
C_0 = -3.11680E-07;

%% Read in & process metadata from log file

% Read in metadata
[num_bottle txt_bottle raw_bottle] = xlsread(AutoBOD_LogFile,'Bottle_metadata');
[num_deploy txt_deploy raw_deploy] = xlsread(AutoBOD_LogFile,'Deployment_metadata');

% Flow fields into cell array AutoBOD_deploy_metadata

% Easy fields first
AutoBOD_deploy_metadata.Deployment_ID = cellstr(raw_deploy(10:end,1));
AutoBOD_deploy_metadata.Cruise_ID = cellstr(txt_deploy(10:end,2));
AutoBOD_deploy_metadata.Incu_temp_deg_C = cell2mat(raw_deploy(10:end,5));
AutoBOD_deploy_metadata.Datafile_loc = cellstr(raw_deploy(10:end,6));
AutoBOD_deploy_metadata.Presens_ores = cell2mat(raw_deploy(10:end,8));
AutoBOD_deploy_metadata.Presens_oxyu = cell2mat(raw_deploy(10:end,7));
AutoBOD_deploy_metadata.Baro_press_hPa = cell2mat(raw_deploy(10:end,9));
AutoBOD_deploy_metadata.Computerclock_tz = cellstr(raw_deploy(10:end,12));
AutoBOD_bottle_metadata.Deployment_ID = cellstr(raw_bottle(6:end,1));
AutoBOD_bottle_metadata.Bottle_ID = cell2mat(raw_bottle(6:end,2));
AutoBOD_bottle_metadata.Sample_ID = cellstr(raw_bottle(6:end,3));
AutoBOD_bottle_metadata.CTD_ID = cellstr(txt_bottle(6:end,4));
AutoBOD_bottle_metadata.Sample_replicate = cell2mat(raw_bottle(6:end,5));
AutoBOD_bottle_metadata.Depth_m = cell2mat(raw_bottle(6:end,6));
AutoBOD_bottle_metadata.Salinity = cell2mat(raw_bottle(6:end,7));
AutoBOD_bottle_metadata.Lat = cell2mat(raw_bottle(6:end,8));
AutoBOD_bottle_metadata.Long = cell2mat(raw_bottle(6:end,9));

% Get number of total deployments in log
NumDeploys = length(AutoBOD_deploy_metadata.Deployment_ID);

% Now, the fields requiring a bit more manipulation

% Time zone
AutoBOD_deploy_metadata.Incu_timezone = regexp(raw_deploy(10:end,3),'.{6}$','match');
AutoBOD_deploy_metadata.Incu_timezone = ...
    reshape([AutoBOD_deploy_metadata.Incu_timezone{:}],NumDeploys,1);

% Incubation start/endtimes
AutoBOD_deploy_metadata.Incu_starttime_UTC = datetime(raw_deploy(10:end,3),...
    'InputFormat','uuuu-MM-dd''T''HH:mm:ssXXX','TimeZone','UTC');
AutoBOD_deploy_metadata.Incu_endtime_UTC = datetime(raw_deploy(10:end,4),...
    'InputFormat','uuuu-MM-dd''T''HH:mm:ssXXX','TimeZone','UTC');

% Start/endtimes to be used for calculation of respiration rates (if given)

for i=1:length(AutoBOD_deploy_metadata.Deployment_ID)
    
    if strcmp([txt_deploy(9+i,10)],'')
        
        % No specific calculation timepoints were specified, just use
        % deployment start/end points
        
        AutoBOD_deploy_metadata.Calc_starttime_UTC(i) = AutoBOD_deploy_metadata.Incu_starttime_UTC(i);
        AutoBOD_deploy_metadata.Calc_endtime_UTC(i) = AutoBOD_deploy_metadata.Incu_endtime_UTC(i);
        
    else
        
        % The user specified some calculation start/end points
        
        AutoBOD_deploy_metadata.Calc_starttime_UTC(i) = datetime(raw_deploy(9+i,10),...
            'InputFormat','uuuu-MM-dd''T''HH:mm:ssXXX','TimeZone','UTC');
        AutoBOD_deploy_metadata.Calc_endtime_UTC(i) = datetime(raw_deploy(9+i,11),...
            'InputFormat','uuuu-MM-dd''T''HH:mm:ssXXX','TimeZone','UTC');
        
    end
    
end

clear i;

% Create some other timestamps in local time (relative to where deployment
% took place)

% Not creating these as datetime objects, since it appears the entire
% object has to have the same time zone associated with it; this not
% acceptable for our purposes since many of the deployments were conducted
% in different time zones

% Deployment times

for i=1:length(AutoBOD_deploy_metadata.Deployment_ID)
    AutoBOD_deploy_metadata.Incu_starttime_local{i} = ...
        datestr(datetime(AutoBOD_deploy_metadata.Incu_starttime_UTC(i),'TimeZone',...
        AutoBOD_deploy_metadata.Incu_timezone{i}));
    AutoBOD_deploy_metadata.Incu_endtime_local{i} = ...
        datestr(datetime(AutoBOD_deploy_metadata.Incu_endtime_UTC(i),'TimeZone',...
        AutoBOD_deploy_metadata.Incu_timezone{i}));
end

clear i

% Times to be used for calculation of rates

for i=1:length(AutoBOD_deploy_metadata.Deployment_ID)
    AutoBOD_deploy_metadata.Calc_starttime_local{i} = ...
        datestr(datetime(AutoBOD_deploy_metadata.Calc_starttime_UTC(i),'TimeZone',...
        AutoBOD_deploy_metadata.Incu_timezone{i}));
    AutoBOD_deploy_metadata.Calc_endtime_local{i} = ...
        datestr(datetime(AutoBOD_deploy_metadata.Calc_endtime_UTC(i),'TimeZone',...
        AutoBOD_deploy_metadata.Incu_timezone{i}));
end

clear i

% Calculate deployment durations
AutoBOD_deploy_metadata.Deploy_duration_hours = ...
    AutoBOD_deploy_metadata.Incu_endtime_UTC - ...
    AutoBOD_deploy_metadata.Incu_starttime_UTC;

% Calculate durations of data range over which rates are to be calculated
AutoBOD_deploy_metadata.Calcrange_duration_hours = ...
    AutoBOD_deploy_metadata.Calc_endtime_UTC - ...
    AutoBOD_deploy_metadata.Calc_starttime_UTC;

% Create vector of deployments to be analyzed
if (exist('Desired_Deployments') & ~isempty('Desired_Deployments'))
    Deploy_queue = cellstr(Desired_Deployments);
else
    Deploy_queue = unique(AutoBOD_deploy_metadata.Deployment_ID);
end

% Finally, preallocate some elements of a structure AutoBOD_rate_results
% into which we will flow our rate calculations

AutoBOD_rate_results.Deployment_ID = num2cell(NaN(length(AutoBOD_bottle_metadata.Deployment_ID),1));
AutoBOD_rate_results.Bottle_ID = NaN(length(AutoBOD_bottle_metadata.Deployment_ID),1);
AutoBOD_rate_results.Sample_ID = num2cell(NaN(length(AutoBOD_bottle_metadata.Deployment_ID),1));
AutoBOD_rate_results.Sample_replicate = NaN(length(AutoBOD_bottle_metadata.Deployment_ID),1);
AutoBOD_rate_results.dO2dt_umol_L_d = NaN(length(AutoBOD_bottle_metadata.Deployment_ID),1);
AutoBOD_rate_results.dO2dt_umol_L_d_sigma = NaN(length(AutoBOD_bottle_metadata.Deployment_ID),1);

%% Analysis, for each desired deployment in Deployment_queue

for i=1:length(Deploy_queue)
    
%     %% Pull out bottle metadata for this deployment
%     % Doesn't look like we need this anymore (9/10/15)
%     
%     % Store in AutoBOD_bottle_metadata_thisdeploy
%     
%     AutoBOD_bottle_metadata_thisdeploy.Deployment_ID = ...
%         AutoBOD_bottle_metadata.Deployment_ID(find(strcmp([...
%         AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i})));
%     AutoBOD_bottle_metadata_thisdeploy.Bottle_ID = ...
%         AutoBOD_bottle_metadata.Bottle_ID(find(strcmp([...
%         AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i})));
%     AutoBOD_bottle_metadata_thisdeploy.Sample_ID = ...
%         AutoBOD_bottle_metadata.Sample_ID(find(strcmp([...
%         AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i})));
%     AutoBOD_bottle_metadata_thisdeploy.CTD_ID = ...
%         AutoBOD_bottle_metadata.CTD_ID(find(strcmp([...
%         AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i})));
%     AutoBOD_bottle_metadata_thisdeploy.Sample_replicate = ...
%         AutoBOD_bottle_metadata.Sample_replicate(find(strcmp([...
%         AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i})));
%     AutoBOD_bottle_metadata_thisdeploy.Depth_m = ...
%         AutoBOD_bottle_metadata.Depth_m(find(strcmp([...
%         AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i})));
%     AutoBOD_bottle_metadata_thisdeploy.Salinity = ...
%         AutoBOD_bottle_metadata.Salinity(find(strcmp([...
%         AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i})));
%     AutoBOD_bottle_metadata_thisdeploy.Lat = ...
%         AutoBOD_bottle_metadata.Lat(find(strcmp([...
%         AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i})));
%     AutoBOD_bottle_metadata_thisdeploy.Long = ...
%         AutoBOD_bottle_metadata.Long(find(strcmp([...
%         AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i})));
    
    %% Read-in and parsing of data file
    
    % First, open the datafile for this deployment and determine number of
    % header lines, i.e., initial lines of numbers in the data file that aren't
    % "good" AutoBOD data
    
    HeadCount = 0; % Create a variable to keep track of number of header lines
    
    % Get index of this deployment in metadata array
    Ind_ThisDeploy = find(strcmp([...
        AutoBOD_deploy_metadata.Deployment_ID'],Deploy_queue{i}));
    
    % Open the datafile for this deployment
    FileID = fopen(AutoBOD_deploy_metadata.Datafile_loc{Ind_ThisDeploy});

    while (~feof(FileID))
        TLine = fgetl(FileID); % Get one line of data at a time and evaluate
        if ~isempty(regexp(TLine,'\d{5}\s(\d{6}|[-]\d{5})\s\d{2}[.]\d{2}\s\d{3}[.]\d{2}\s\d\s\d\s\d{5}\s\d{3}\s\d\s\d{2}/\d{2}/\d{2}\s\d{2}:\d{2}:\d{2}'))
            % We've hit good data matching the desired pattern, so we can stop
            break;
        end
        HeadCount = HeadCount + 1; % Otherwise, increment our counter
    end
    
    fclose(FileID);
    
    % Now that we know number of lines to ignore, can use textscan to read
    % everything in
   
    % Open the datafile for this deployment (again)
    FileID = fopen(AutoBOD_deploy_metadata.Datafile_loc{Ind_ThisDeploy});
    
    % Read in data using textscan; by its nature, will stop at either EOF
    % or when it hits a line of stuff that doesn't match the pattern
    
    AutoBOD_rawdata = textscan(FileID,'%f %f %f %f %d %d %f %f %f %s %s','HeaderLines',HeadCount);

    fclose(FileID);
    
    % Now, deal with each field
        
    % ****** Create timestamps ******
    
    % Have to first account for computer clock, then convert to whatever
    % timezone the deployment actually occured in
    
    AutoBOD_data.Timestamp_local = datetime(year(datenum(AutoBOD_rawdata{10})),...
        month(datenum(AutoBOD_rawdata{10})),...
        day(datenum(AutoBOD_rawdata{10})),...
        hour(datenum(AutoBOD_rawdata{11})),...
        minute(datenum(AutoBOD_rawdata{11})),...
        second(datenum(AutoBOD_rawdata{11})),'TimeZone',...
        AutoBOD_deploy_metadata.Computerclock_tz{Ind_ThisDeploy});
    
    AutoBOD_data.Timestamp_local = datetime(AutoBOD_data.Timestamp_local,'TimeZone',...
        AutoBOD_deploy_metadata.Incu_timezone{Ind_ThisDeploy});

    % ****** Temperature ****** 
    
    AutoBOD_data.Temp_deg_C = AutoBOD_rawdata{3};
    
    % ****** DO ******
    
    % Read in raw, salinity uncorrected value
    
    AutoBOD_data.DO_uncorr = AutoBOD_rawdata{4};
    
    % Now, convert to umol/L with decimal in right place; this depends on
    % the units/format of the values recorded during the deployment
    
    % Note that we will correct for salinity later, once we figure out
    % what data segment came from which bottle
    
    % ***** This section works ok for data in % sat with ores = 2,
    % but user should revisit this section after making any changes to OEM
    % settings such that data input format is different
    
    if AutoBOD_deploy_metadata.Presens_oxyu(Ind_ThisDeploy)==0
        
        % Code oxyu = 0 indicates oxygen was obtained in % sat
        
        % First, put decimal in correct place
        
        AutoBOD_data.DO_uncorr_dec_adj = AutoBOD_data.DO_uncorr*...
            (10^(AutoBOD_deploy_metadata.Presens_ores(Ind_ThisDeploy)-1));
        
        % Now, convert from % sat to uM using same Bunsen coefficients etc.
        % as PreSens
        
        AtmP = AutoBOD_deploy_metadata.Baro_press_hPa(Ind_ThisDeploy); % Retrieve atmos. press. in hPa
        
        AutoBOD_data.DO_uM_O2_uncorr = ((AtmP-exp(52.57-6690.9./(273.15+AutoBOD_data.Temp_deg_C)-4.681.*log(273.15+AutoBOD_data.Temp_deg_C)))/1013)...
            .*AutoBOD_data.DO_uncorr_dec_adj./100*0.2095.*...
            (48.998-1.335.*AutoBOD_data.Temp_deg_C+0.02755.*...
            power(AutoBOD_data.Temp_deg_C,2)-...
            0.000322.*power(AutoBOD_data.Temp_deg_C,3)+...
            0.000001598.*power(AutoBOD_data.Temp_deg_C,4))*32/22.414*31.25;
        
    elseif AutoBOD_deploy_metadata.Presens_oxyu(Ind_ThisDeploy)==5
            
        % Code oxyu = 5 indicates oxygen was obtained in micromoles/L
        
        % First, put decimal in correct place
        
        AutoBOD_data.DO_uncorr_dec_adj = AutoBOD_data.DO_uncorr*...
            (10^(AutoBOD_deploy_metadata.Presens_ores(Ind_ThisDeploy)-1));
        
        % Shouldn't require any conversion, so take data as it is now
        
         AutoBOD_data.DO_uM_O2_uncorr = AutoBOD_data.DO_uncorr_dec_adj;      
        
    else
        
        % Future code to convert data obtained in other formats should go
        % here
        
        AutoBOD_data.DO_uM_O2_uncorr = AutoBOD_data.DO_uncorr;
        
    end
    
    % ****** Error ****** 
    
    AutoBOD_data.Err = AutoBOD_rawdata{5};
  
    % ****** Bottle lock ****** 
    
    AutoBOD_data.Bot_lock = AutoBOD_rawdata{6};    
            
    % ****** Home detect code ****** 
    
    AutoBOD_data.HD = AutoBOD_rawdata{9};
    
    %% Association of each line of DO data in this dataset with a particular bottle number
    
    % Cumbersome, but will accomplish by cycling through each line of data
    % to see if it matches the desired pattern
    
    % Note that we assume the first bottle with data was bottle #1
    
    % Since the data output from the instrument isn't perfect, will have to
    % use multiple criteria to try and hone in on the "good" segments of
    % data
    
    % Often, this will result in some conservative "clipping" of a few
    % lines of good data on either end of a segment
    
    % First, initialize vectors to hold to the bottle numbers we're about
    % to assign, and a sequential segment ID for each segment of good data
    
    AutoBOD_data.Bot_ID = NaN(length(AutoBOD_data.DO_uM_O2_uncorr),1);
    AutoBOD_data.Segment = NaN(length(AutoBOD_data.DO_uM_O2_uncorr),1);

    % Also, create some counters to keep track of what kind of data
    % we're dealing with
    
    Badcount = 0;
    This_bottle = 1;
    This_segment = 1;

    % Now, cycle through the data
    
    for j=1:length(AutoBOD_data.DO_uM_O2_uncorr)
        
        if AutoBOD_data.DO_uM_O2_uncorr(j) > 800
            
            % As a first check, make sure we don't have a biologically
            % implausible value (the instrument spits out 999.999 as a
            % value in the absence of a good reading)
            
            Badcount = Badcount + 1;
            
        else
            
            if (AutoBOD_data.Err(j)==0 && AutoBOD_data.Bot_lock(j)==1)
                
                % We probably have a line of good data, but have to make
                % sure it isn't an isolated, anomalous reading
                % (these exist in some of the data)
                
                % Define a window in which we'll consider the quality of
                % the data
                
                Window_low = j-Max_RECD;
                
                if j<length(AutoBOD_data.DO_uM_O2_uncorr)-Max_RECD
                    Window_hi = j+Max_RECD;
                else
                    Window_hi = j+(length(AutoBOD_data.DO_uM_O2_uncorr)-j);
                end
                
                if (sum(AutoBOD_data.Err(Window_low:Window_hi))<=Max_RECD && ...
                        sum(AutoBOD_data.Bot_lock(Window_low:Window_hi))>=Max_RECD*1.5)
                    
                    % This line is almost certainly one good line of data within a
                    % series of lines of good data; the cutoffs above
                    % were tuned based on the quantity of "sporadic" bad readings
                    % observed in a number of datasets, but you can still
                    % change them if you want
                
                    AutoBOD_data.Bot_ID(j) = This_bottle;
                    AutoBOD_data.Segment(j) = This_segment;
                    Badcount = 0;
                    
                else
                    
                    Badcount = Badcount + 1;
                    
                end
                    
            else % It's not an "error-free" line of data, so we have to do some investigating
                
                if (j<=Max_RECD || sum(AutoBOD_data.Err(j-Max_RECD:j))==0)
                    
                    % We can be certain we haven't had an error code in the past few
                    % reads
                     
                    if (AutoBOD_data.Err(j)==0 && AutoBOD_data.Bot_lock(j)==0 && ...
                            Badcount<=Max_RECD && j>=3)
                        
                        % We'll give this read the benefit of the doubt, but note it as a
                        % potentially bad line (so long as we've had a recent good
                        % read)
                        
                        AutoBOD_data.Bot_ID(j) = This_bottle;
                        Badcount = Badcount + 1;
                        
                    elseif (AutoBOD_data.Err(j)==0 && AutoBOD_data.Bot_lock(j)==0 && ...
                            Badcount>Max_RECD)
                        
                        % We're likely into bad data
                        
                        Badcount = Badcount + 1;
                        
                    elseif AutoBOD_data.Err(j)==4
                        
                        % We definitely have bad data
                        
                        Badcount = Badcount + 1;
                        
                    end
                    
                else
                    
                    Badcount = Badcount + 1;
                    
                end
                
            end
                        
        end
        
        % Before moving onto the next line of data, we have to determine whether
        % we're still on the same bottle number
        
        if (Badcount==Max_RECD+1 && j>Max_RECD+1)
            
            % Safe to conclude we're at the end of the data for this
            % particular bottle, let's get ready for the next one
            
            % Increment bottle ID
            
            if This_bottle<=(Num_bots-1)
                
                This_bottle = This_bottle+1;
                
            elseif This_bottle==Num_bots
                
                This_bottle = 1;
                
            end
            
            % Increment segment ID
            
            This_segment = This_segment+1;
             
        end
            
    end
    
    clear j;
    
    %% Calculation of data segment averages, with bottle ID and timestamp
    
    % Includes correction of data for salinity, since we now know which
    % segments of data go with which bottles
    
    % Before beginning, preallocate vectors in a structure AutoBOD_segmeans
    
    NumSegs = length(unique(AutoBOD_data.Segment(~isnan(AutoBOD_data.Segment))));
    
    AutoBOD_segmeans.Segment = NaN(NumSegs,1);
    AutoBOD_segmeans.Bot_ID = NaN(NumSegs,1);
    AutoBOD_segmeans.Temp_deg_C = NaN(NumSegs,1);
    AutoBOD_segmeans.DO_uM_O2_mean = NaN(NumSegs,1);
    AutoBOD_segmeans.DO_uM_O2_sigma = NaN(NumSegs,1);

    % Timestamp object can't be preallocated in same way, so we'll use
    % existing time data for this dataset as a start
    
    AutoBOD_segmeans.Timestamp_local = AutoBOD_data.Timestamp_local(1:NumSegs);

    Seg_count = 1; % Define a counter so we know where to insert data into AutoBOD_segmeans
    
    for k=1:Num_bots
        
        % First, get the unique segment IDs associated with this bottle ID
        
        Segs_this_bottle = unique(AutoBOD_data.Segment(AutoBOD_data.Bot_ID==k));
        Segs_this_bottle = Segs_this_bottle(~isnan(Segs_this_bottle));
    
        % Also, capture the salinity recorded for this bottle in the
        % metadata log
        
        Sal_this_bottle = AutoBOD_bottle_metadata.Salinity(find(...
                (strcmp([AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i}))' &...
                AutoBOD_bottle_metadata.Bottle_ID==k));
        
        % Now, cycle through segments for this bottle and calculate mean
        % oxygen concentration; record this and other relevant data in
        % AutoBOD_segmeans
        
        for l=1:length(Segs_this_bottle)
            
            % Get indices to data for this bottle and segment
            
            Ind_segdata = find(AutoBOD_data.Segment==Segs_this_bottle(l) &...
                AutoBOD_data.Bot_ID==k);
            
            % ****** Segment number ******
                        
            AutoBOD_segmeans.Segment(Seg_count) = Segs_this_bottle(l);
            
            % ****** Bottle ID ******
            
            AutoBOD_segmeans.Bot_ID(Seg_count) = k;
            
            % ****** Temp ******
            
            AutoBOD_segmeans.Temp_deg_C(Seg_count) = ...
                mean(AutoBOD_data.Temp_deg_C(Ind_segdata));

            % ****** Timestamp ******
            
            % Using timepoint of the central datapoint in this segment
            
            AutoBOD_segmeans.Timestamp_local(Seg_count) = ...
                AutoBOD_data.Timestamp_local(round(median(Ind_segdata)));
            
            % ****** Salinity-corrected DO, in uM ******
            
            % Correct for salinity, now that we know what bottle number
            % goes with what data, per:
            
            % Garcia HE, Gordon LI (1992), Oxygen solubility in seawater: Better
            % fitting equations, Limnol. Oceanogr. 37:1307-1312
            
            % with coefficients as defined above
            
            % ***** Critical: Assumes a salinity setting of "0" was used when data
            % were collected; no need to correct for pressure since we're
            % conducting these incubations in lab
            
            % Calculate a scaled temperature for each reading
            
            Scaled_temps = log((298.15-AutoBOD_data.Temp_deg_C(Ind_segdata))./...
                (273.15+AutoBOD_data.Temp_deg_C(Ind_segdata)));
            
            % Calculate a salinity compensation factor
            
            Sal_comp_factors = exp(Sal_this_bottle.*(...
                B_0+...
                B_1.*Scaled_temps+...
                B_2.*(Scaled_temps.^2)+...
                B_3.*(Scaled_temps.^3))+...
                C_0.*(Sal_this_bottle^2));
            
            % Compute salinity-corrected DO concentrations in umol/L
            
            DO_uM_O2_salcor = AutoBOD_data.DO_uM_O2_uncorr(Ind_segdata).*Sal_comp_factors;
            
            % Calculate mean and sd and store
            
            AutoBOD_segmeans.DO_uM_O2_mean(Seg_count) = mean(DO_uM_O2_salcor);
            AutoBOD_segmeans.DO_uM_O2_sigma(Seg_count) = std(DO_uM_O2_salcor);
            
            % Lastly, advance our counter
            
            Seg_count = Seg_count+1;
            
        end
        
    end
    
    clear k l;
    
    %% Plot data for this deployment
    
    figure;

    Plot_colors=hsv(Num_bots);
    Legend_text = num2cell(NaN(Num_bots,1));

    hold on;
    for n=1:Num_bots
        
        % Plot data
        
        plot(AutoBOD_segmeans.Timestamp_local(AutoBOD_segmeans.Bot_ID==n),...
            AutoBOD_segmeans.DO_uM_O2_mean(AutoBOD_segmeans.Bot_ID==n),...
            'Color',Plot_colors(n,:),'LineStyle','-','Marker','none');
        Legend_text{n} = ['Bottle ' num2str(n) ', ' strrep(AutoBOD_bottle_metadata.Sample_ID{find(...
            (strcmp([AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i}))' &...
            AutoBOD_bottle_metadata.Bottle_ID==n)},'_','\_')];
    end
    
    % Superimpose lines showing segment used for rate calculation
    
    % Retrieve start and end times to be used for calculation in next
    % section
    
    t1 = AutoBOD_deploy_metadata.Calc_starttime_local(Ind_ThisDeploy);
    t2 = AutoBOD_deploy_metadata.Calc_endtime_local(Ind_ThisDeploy);
    tz = AutoBOD_deploy_metadata.Incu_timezone{Ind_ThisDeploy};
        
    plot([datetime(t1,'TimeZone',tz) datetime(t1,'TimeZone',tz)],ylim,'k--')
    plot([datetime(t2,'TimeZone',tz) datetime(t2,'TimeZone',tz)],ylim,'k--')
    
    title(['AutoBOD data for deployment: ' strrep(Deploy_queue{i},'_','\_')]);
    xlabel('Time');
    ylabel('Dissolved oxygen (\mumol/L)');
    legend(Legend_text);

    hold off;
    
    clear n;

    %% Calculation of rates of DO change in each bottle, using linear regression
    
    % Requires linfit.m
    
    for m=1:Num_bots
        
        % Extract needed data for this bottle
        
        % Retrieve start and end times to be used for calculation
        
        t1 = AutoBOD_deploy_metadata.Calc_starttime_local(Ind_ThisDeploy);
        t2 = AutoBOD_deploy_metadata.Calc_endtime_local(Ind_ThisDeploy);
        
        x = juliandate(AutoBOD_segmeans.Timestamp_local(AutoBOD_segmeans.Bot_ID==m &...
            AutoBOD_segmeans.Timestamp_local>=t1 &...
            AutoBOD_segmeans.Timestamp_local<=t2));
        y = AutoBOD_segmeans.DO_uM_O2_mean(AutoBOD_segmeans.Bot_ID==m &...
            AutoBOD_segmeans.Timestamp_local>=t1 &...
            AutoBOD_segmeans.Timestamp_local<=t2);
        sy = AutoBOD_segmeans.DO_uM_O2_sigma(AutoBOD_segmeans.Bot_ID==m &...
            AutoBOD_segmeans.Timestamp_local>=t1 &...
            AutoBOD_segmeans.Timestamp_local<=t2);
        
        [a sa cov r]=linfit(x,y,sy);
        
        % Record rate and error in AutoBOD_rate_results
        
        AutoBOD_rate_results.dO2dt_umol_L_d(find(...
            (strcmp([AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i}))' &...
            AutoBOD_bottle_metadata.Bottle_ID==m)) = a(2);
        
         AutoBOD_rate_results.dO2dt_umol_L_d_sigma(find(...
            (strcmp([AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i}))' &...
            AutoBOD_bottle_metadata.Bottle_ID==m)) = sa(2);
        
        % Record metadata
        
        AutoBOD_rate_results.Deployment_ID(find(...
            (strcmp([AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i}))' &...
            AutoBOD_bottle_metadata.Bottle_ID==m)) = cellstr(Deploy_queue{i});
        
        AutoBOD_rate_results.Bottle_ID(find(...
            (strcmp([AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i}))' &...
            AutoBOD_bottle_metadata.Bottle_ID==m)) = m;
        
        AutoBOD_rate_results.Sample_ID(find(...
            (strcmp([AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i}))' &...
            AutoBOD_bottle_metadata.Bottle_ID==m)) = cellstr(AutoBOD_bottle_metadata.Sample_ID(find(...
            (strcmp([AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i}))' &...
            AutoBOD_bottle_metadata.Bottle_ID==m)));
        
        AutoBOD_rate_results.Sample_replicate(find(...
            (strcmp([AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i}))' &...
            AutoBOD_bottle_metadata.Bottle_ID==m)) = AutoBOD_bottle_metadata.Sample_replicate(find(...
            (strcmp([AutoBOD_bottle_metadata.Deployment_ID'],Deploy_queue{i}))' &...
            AutoBOD_bottle_metadata.Bottle_ID==m));
        
        
    end
        
end

%% Finally, can group replicate bottle data by treatment & calculate overall rates

% Preallocate elements of destination structure 

% Get unique deployment-treatment combinations that exist in the metadata
% log

TreatsTable = table(AutoBOD_rate_results.Deployment_ID,AutoBOD_rate_results.Sample_ID);
UniqueTreats = unique(TreatsTable);

NumCombs = height(UniqueTreats);

AutoBOD_treat_results.Deployment_ID = UniqueTreats{:,1};
AutoBOD_treat_results.Sample_ID = UniqueTreats{:,2};
AutoBOD_treat_results.NumReps = NaN(NumCombs,1);
AutoBOD_treat_results.dO2dt_umol_L_d = NaN(NumCombs,1);
AutoBOD_treat_results.dO2dt_umol_L_d_sigma = NaN(NumCombs,1);

for i=1:NumCombs
    
    AutoBOD_treat_results.NumReps(i) = ...
        length(TreatsTable{strcmp(TreatsTable{:,1},UniqueTreats{i,1}) &...
        strcmp(TreatsTable{:,2},UniqueTreats{i,2}),:});
    
    AutoBOD_treat_results.dO2dt_umol_L_d(i) = ...
        mean(AutoBOD_rate_results.dO2dt_umol_L_d(strcmp(AutoBOD_rate_results.Deployment_ID,AutoBOD_treat_results.Deployment_ID(i)) &...
        strcmp(AutoBOD_rate_results.Sample_ID,AutoBOD_treat_results.Sample_ID(i))));
    
    AutoBOD_treat_results.dO2dt_umol_L_d_sigma(i) = sqrt(sum(AutoBOD_rate_results.dO2dt_umol_L_d_sigma(strcmp(AutoBOD_rate_results.Deployment_ID,AutoBOD_treat_results.Deployment_ID(i)) &...
        strcmp(AutoBOD_rate_results.Sample_ID,AutoBOD_treat_results.Sample_ID(i)))))/AutoBOD_treat_results.NumReps(i);
    
end
    
    
    