%% PHORCYS_data_analysis_Model_B.m

% Owner/author: James Collins, MIT-WHOI Joint Program, Woods Hole
% Oceanographic Institution, james.r.collins@aya.yale.edu

% Process data from the PHORCYS Model "B" (magnetic drive repeatedly opens
% and closes chambers) and calculate community metabolic rates

% The Model "B" instrument uses Aanderaa Aqua 4531 optodes

% Created 11/11/16 by JRC under MATLAB R2015a; maintained on GitHub at
% https://github.com/jamesrco/DO_Instruments under a GNU General Public
% License v. 3

% All time values are in local time

%% -------------------------------------------------------------------------

% clear & close to prep workspace

clear all;

close all;

% set constants for salinity compensation, see Aanderaa 4330/4835 TD269 manual p. 32

B_0 = -6.24097e-3;
B_1 = -6.93498e-3;
B_2 = -6.90358e-3;
B_3 = -4.29155e-3;
C_0 = -3.11680e-7;


% -------------------------------------------------------------------------
% load PHORCYS data from model "B" instrument

% ** these data have a different format **

% ** assumes user has manually concatenated 4-line data strings at each
% timepoint into a single line **
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% read in data for Iselin dock (WHOI) deployments in November 2016
% -------------------------------------------------------------------------

% the data exist in separate files for each time interval following a
% chamber opening; openings were at 0600 and 1700 local time, to more or
% less coincide with sunrise and sunset

% define a list of dates and times (corresponding to filenames)

Iselin_201611

% open the datafile for this deployment

Iselin_20161107_1111 = csvread('/Users/jrcollins/Code/DO_Instruments/PHORCYS/data/processed/Iselin_WHOI_2016_11/PHORCYS_110716_110816.csv');
    
KN207_3_PS1_20120617_Timestamp_DO = KN207_3_PS1_20120617(:,1);
KN207_3_PS1_20120617_Timestamp_DO = x2mdate(KN207_3_PS1_20120617_Timestamp_DO,1);
KN207_3_PS1_20120617_timefrac = KN207_3_PS1_20120617(:,2);
KN207_3_PS1_20120617_DO_1_dark_uM_raw = KN207_3_PS1_20120617(:,3);
KN207_3_PS1_20120617_DO_1_dark_sat = KN207_3_PS1_20120617(:,5);
KN207_3_PS1_20120617_DO_1_dark_temp_degC = KN207_3_PS1_20120617(:,6);
KN207_3_PS1_20120617_DO_2_light_uM_raw = KN207_3_PS1_20120617(:,7);
KN207_3_PS1_20120617_DO_2_light_sat = KN207_3_PS1_20120617(:,8);
KN207_3_PS1_20120617_DO_2_light_temp_degC = KN207_3_PS1_20120617(:,9);

% adjust DO values for salinity, see Aanderaa manual p. 32

KN207_3_PS1_20120617_DO_1_dark_uM = zeros(length(KN207_3_PS1_20120617_DO_1_dark_uM_raw),1);

for i=1:length(KN207_3_PS1_20120617_DO_1_dark_uM_raw)
    closest_sal = find(KN207_3_mettime>=KN207_3_PS1_20120617_Timestamp_DO(i),1,'first');
    sal = KN207_3_sal(closest_sal);
    T_s = log((298.15-KN207_3_PS1_20120617_DO_1_dark_temp_degC)./(273.15+KN207_3_PS1_20120617_DO_1_dark_temp_degC));
    KN207_3_PS1_20120617_DO_1_dark_uM=KN207_3_PS1_20120617_DO_1_dark_uM_raw.*exp(sal*(B_0+B_1*T_s+B_2*(T_s).^2+B_3*(T_s).^3)+C_0*sal.^2);
end

KN207_3_PS1_20120617_DO_2_light_uM = zeros(length(KN207_3_PS1_20120617_DO_2_light_uM_raw),1);

for i=1:length(KN207_3_PS1_20120617_DO_2_light_uM_raw)
    closest_sal = find(KN207_3_mettime>=KN207_3_PS1_20120617_Timestamp_DO(i),1,'first');
    sal = KN207_3_sal(closest_sal);
    T_s = log((298.15-KN207_3_PS1_20120617_DO_2_light_temp_degC)./(273.15+KN207_3_PS1_20120617_DO_2_light_temp_degC));
    KN207_3_PS1_20120617_DO_2_light_uM=KN207_3_PS1_20120617_DO_2_light_uM_raw.*exp(sal*(B_0+B_1*T_s+B_2*(T_s).^2+B_3*(T_s).^3)+C_0*sal.^2);
end

% load Winkler data

KN207_3_PS1_dark_Wink_time=xlsread('/Users/jrcollins/Code/DO_Instruments/PHORCYS/data/processed/KN207_3/KN207_3_Process_Station_1_PHORCYS_2012_06_17_19.xlsx','Winklers','A8');
KN207_3_PS1_dark_Wink_time = x2mdate(KN207_3_PS1_dark_Wink_time,1);
KN207_3_PS1_dark_Wink_uM=xlsread('/Users/jrcollins/Code/DO_Instruments/PHORCYS/data/processed/KN207_3/KN207_3_Process_Station_1_PHORCYS_2012_06_17_19.xlsx','Winklers','I8');

% normalize light bottle data to calibrated dark bottle data

darklightdiff=KN207_3_PS1_20120617_DO_1_dark_uM-KN207_3_PS1_20120617_DO_2_light_uM;
meandarklightdiff=mean(darklightdiff(96:260,:));

% calibrate dark data to Winklers

Wink_caltime = find(KN207_3_PS1_20120617_Timestamp_DO>=KN207_3_PS1_dark_Wink_time,1,'first');
Wink_diff = KN207_3_PS1_dark_Wink_uM-KN207_3_PS1_20120617_DO_1_dark_uM(Wink_caltime);
Wink_diff_19Jun=Wink_diff;
KN207_3_PS1_20120617_DO_1_dark_uM_cal = KN207_3_PS1_20120617_DO_1_dark_uM+Wink_diff;
KN207_3_PS1_20120617_DO_2_light_uM_cal=(KN207_3_PS1_20120617_DO_2_light_uM+meandarklightdiff)+Wink_diff;

% truncate data set to include only usable data

KN207_3_PS1_20120617_data = [KN207_3_PS1_20120617_Timestamp_DO KN207_3_PS1_20120617_timefrac ...
    KN207_3_PS1_20120617_DO_1_dark_uM_cal KN207_3_PS1_20120617_DO_1_dark_temp_degC ...
    KN207_3_PS1_20120617_DO_2_light_uM_cal KN207_3_PS1_20120617_DO_2_light_temp_degC];

% KN207_3_PS1_20120617_data = KN207_3_PS1_20120617_data(142:1313,:);

KN207_3_PS1_20120617_data = KN207_3_PS1_20120617_data(39:end,:);

KN207_3_PS1_20120617_Timestamp_DO = KN207_3_PS1_20120617_data(:,1);
KN207_3_PS1_20120617_timefrac = KN207_3_PS1_20120617_data(:,2);
KN207_3_PS1_20120617_DO_1_dark_uM = KN207_3_PS1_20120617_data(:,3);
KN207_3_PS1_20120617_DO_1_dark_temp_degC = KN207_3_PS1_20120617_data(:,4);
KN207_3_PS1_20120617_DO_2_light_uM = KN207_3_PS1_20120617_data(:,5);
KN207_3_PS1_20120617_DO_2_light_temp_degC = KN207_3_PS1_20120617_data(:,6);