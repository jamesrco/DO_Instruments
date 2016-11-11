%% PHORCYS_data_analysis_Model_A.m

% Process data from the PHORCYS Model "A" ("one and done" burn wire closure
% system) instrument and calculate community metabolic rates

% All of the deployments processed by this script were done on cruises
% aboard the R/V Knorr; met & ship data are from shipboard instruments
% aboard the Knorr

% The Model "A" instrument collected data from Aanderaa 4330F optodes
% w/default salinity set to 0

% Original script created 7/23/13 by JRC; modf'd 9/4/14, 10/16/15 JRC

% Version history subsequent to 10/16/15 recorded on GitHub at 
% https://github.com/jamesrco/DO_Instruments

% Reuse: GNU General Public License v. 3

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
% read in met data & clean up PAR series
% -------------------------------------------------------------------------

% read in KN207-1 met data

KN207_1_metdata = xlsread('/Users/jrcollins/Code/DO_Instruments/PHORCYS/data/processed/KN207_1/KN207_1_met_data.xlsx');
KN207_1_mettime = KN207_1_metdata(:,3);
KN207_1_mettime = x2mdate(KN207_1_mettime,1);
KN207_1_SPAR_uE_cm2_s = KN207_1_metdata(:,4);
KN207_1_sal = KN207_1_metdata(:,5);

% clean up PAR data

KN207_1_SPAR_uE_cm2_s(KN207_1_SPAR_uE_cm2_s<=0) = 0;


% read in KN207-3 met data

KN207_3_metdata = xlsread('/Users/jrcollins/Code/DO_Instruments/PHORCYS/data/processed/KN207_3/KN207_03_alongtrack.xlsx');
KN207_3_mettime = KN207_3_metdata(:,7);
KN207_3_mettime = x2mdate(KN207_3_mettime,1);
KN207_3_SPAR_uE_cm2_s = KN207_3_metdata(:,19);
KN207_3_sal = KN207_3_metdata(:,12);

% clean up PAR data

KN207_3_SPAR_uE_cm2_s(KN207_3_SPAR_uE_cm2_s<=0) = 0;

% -------------------------------------------------------------------------
% load in data for BLATZ II (KN207-1) station QL-1, 24-27 Apr 2012
% -------------------------------------------------------------------------

% rad data not loaded -- instrument was deployed too deeply to see any
% signal in clear bottle

KN207_1_QL1_20120424 = xlsread('/Users/jrcollins/Code/DO_Instruments/PHORCYS/data/processed/KN207_1/KN207_1_QL1_PHORCYS_Drifter_2012_04_24_27.xlsx','Clean data');

KN207_1_QL1_20120424_Timestamp_DO = KN207_1_QL1_20120424(:,6);
KN207_1_QL1_20120424_Timestamp_DO = x2mdate(KN207_1_QL1_20120424_Timestamp_DO,1);
KN207_1_QL1_20120424_timefrac = KN207_1_QL1_20120424(:,7);
KN207_1_QL1_20120424_DO_1_dark_uM_raw = KN207_1_QL1_20120424(:,8);
KN207_1_QL1_20120424_DO_1_dark_temp_degC = KN207_1_QL1_20120424(:,10);

% adjust DO values for salinity, see Aanderaa 4330/4835 TD269 manual p. 32

KN207_1_QL1_20120424_DO_1_dark_uM = zeros(length(KN207_1_QL1_20120424_DO_1_dark_uM_raw),1);

for i=1:length(KN207_1_QL1_20120424_DO_1_dark_uM_raw)
    closest_sal = find(KN207_1_mettime>=KN207_1_QL1_20120424_Timestamp_DO(i),1,'first');
    sal = KN207_1_sal(closest_sal);
    T_s = log((298.15-KN207_1_QL1_20120424_DO_1_dark_temp_degC)./(273.15+KN207_1_QL1_20120424_DO_1_dark_temp_degC));
    KN207_1_QL1_20120424_DO_1_dark_uM=KN207_1_QL1_20120424_DO_1_dark_uM_raw.*exp(sal*(B_0+B_1*T_s+B_2*(T_s).^2+B_3*(T_s).^3)+C_0*sal.^2);
end

% load Winkler data

KN207_1_QL1_dark_Wink_time=xlsread('/Users/jrcollins/Code/DO_Instruments/PHORCYS/data/processed/KN207_1/KN207data_Winklers_oxygen.xlsx','Usefully distilled!','B9');
KN207_1_QL1_dark_Wink_time=addtodate(KN207_1_QL1_dark_Wink_time,-4,'year');
KN207_1_QL1_dark_Wink_time=addtodate(KN207_1_QL1_dark_Wink_time,-1,'day'); 
KN207_1_QL1_dark_Wink_time=addtodate(KN207_1_QL1_dark_Wink_time,-1,'hour'); % Winkler spreadsheet using 1904 date system
KN207_1_QL1_dark_Wink_time = x2mdate(KN207_1_QL1_dark_Wink_time,1);
KN207_1_QL1_dark_Wink_uM=xlsread('/Users/jrcollins/Code/DO_Instruments/PHORCYS/data/processed/KN207_1/KN207data_Winklers_oxygen.xlsx','Usefully distilled!','G9');

% calibrate to Winklers

Wink_caltime = find(KN207_1_QL1_20120424_Timestamp_DO>=KN207_1_QL1_dark_Wink_time,1,'first');
Wink_diff = KN207_1_QL1_dark_Wink_uM-KN207_1_QL1_20120424_DO_1_dark_uM(Wink_caltime);
KN207_1_QL1_20120424_DO_1_dark_uM_cal = KN207_1_QL1_20120424_DO_1_dark_uM+Wink_diff;

% truncate data set to include only usable data

KN207_1_QL1_20120424_data = [KN207_1_QL1_20120424_Timestamp_DO KN207_1_QL1_20120424_timefrac KN207_1_QL1_20120424_DO_1_dark_uM_cal KN207_1_QL1_20120424_DO_1_dark_temp_degC];

% KN207_1_QL1_20120424_data = KN207_1_QL1_20120424_data(249:2188,:);
KN207_1_QL1_20120424_data = KN207_1_QL1_20120424_data(1:end,:);

KN207_1_QL1_20120424_Timestamp_DO = KN207_1_QL1_20120424_data(:,1);
KN207_1_QL1_20120424_timefrac = KN207_1_QL1_20120424_data(:,2);
KN207_1_QL1_20120424_DO_1_dark_uM = KN207_1_QL1_20120424_data(:,3);
KN207_1_QL1_20120424_DO_1_dark_temp_degC = KN207_1_QL1_20120424_data(:,4);

% -------------------------------------------------------------------------
% load in data for BLATZ II (KN207-1) station QL-2, 30 Apr-3 May 2012
% -------------------------------------------------------------------------

% rad data not loaded -- light bottle malfunctioned

KN207_1_QL2_20120430 = xlsread('/Users/jrcollins/Code/DO_Instruments/PHORCYS/data/processed/KN207_1/KN207_1_QL2_PHORCYS_Drifter_2012_04_30_05_03.xlsx','Working data');

KN207_1_QL2_20120430_Timestamp_DO = KN207_1_QL2_20120430(:,2);
KN207_1_QL2_20120430_Timestamp_DO = x2mdate(KN207_1_QL2_20120430_Timestamp_DO,1);
KN207_1_QL2_20120430_timefrac = KN207_1_QL2_20120430(:,3);
KN207_1_QL2_20120430_DO_1_dark_uM_raw = KN207_1_QL2_20120430(:,4);
KN207_1_QL2_20120430_DO_1_dark_sat = KN207_1_QL2_20120430(:,5);
KN207_1_QL2_20120430_DO_1_dark_temp_degC = KN207_1_QL2_20120430(:,6);
% KN207_1_QL2_20120430_DO_2_light_uM_raw = KN207_1_QL2_20120430(:,7);
% KN207_1_QL2_20120430_DO_2_light_sat = KN207_1_QL2_20120430(:,8);
% KN207_1_QL2_20120430_DO_2_light_temp_degC = KN207_1_QL2_20120430(:,9);

% adjust DO values for salinity, see Aanderaa 4330/4835 TD269 manual p. 32

KN207_1_QL2_20120430_DO_1_dark_uM = zeros(length(KN207_1_QL2_20120430_DO_1_dark_uM_raw),1);

for i=1:length(KN207_1_QL2_20120430_DO_1_dark_uM_raw)
    closest_sal = find(KN207_1_mettime>=KN207_1_QL2_20120430_Timestamp_DO(i),1,'first');
    sal = KN207_1_sal(closest_sal);
    T_s = log((298.15-KN207_1_QL2_20120430_DO_1_dark_temp_degC)./(273.15+KN207_1_QL2_20120430_DO_1_dark_temp_degC));
    KN207_1_QL2_20120430_DO_1_dark_uM=KN207_1_QL2_20120430_DO_1_dark_uM_raw.*exp(sal.*(B_0+B_1.*T_s+B_2.*(T_s).^2+B_3.*(T_s).^3)+C_0.*sal.^2);
end

% load Winkler data

KN207_1_QL2_dark_Wink_time=xlsread('/Users/jrcollins/Code/DO_Instruments/PHORCYS/data/processed/KN207_1/KN207oxygen_data_final.xlsx','QL2 data','B7');
KN207_1_QL2_dark_Wink_time = x2mdate(KN207_1_QL2_dark_Wink_time,1);
KN207_1_QL2_dark_Wink_uM=xlsread('/Users/jrcollins/Code/DO_Instruments/PHORCYS/data/processed/KN207_1/KN207oxygen_data_final.xlsx','QL2 data','G7');

% calibrate to Winklers

Wink_diff = KN207_1_QL2_dark_Wink_uM-KN207_1_QL2_20120430_DO_1_dark_uM(2038);
KN207_1_QL2_20120430_DO_1_dark_uM_cal = KN207_1_QL2_20120430_DO_1_dark_uM+Wink_diff;

% truncate data set to include only usable data

KN207_1_QL2_20120430_data = [KN207_1_QL2_20120430_Timestamp_DO KN207_1_QL2_20120430_timefrac KN207_1_QL2_20120430_DO_1_dark_uM_cal KN207_1_QL2_20120430_DO_1_dark_temp_degC];

% KN207_1_QL2_20120430_data = KN207_1_QL2_20120430_data(104:2,:);

KN207_1_QL2_20120430_data = KN207_1_QL2_20120430_data(1:2001,:);

KN207_1_QL2_20120430_Timestamp_DO = KN207_1_QL2_20120430_data(:,1);
KN207_1_QL2_20120430_timefrac = KN207_1_QL2_20120430_data(:,2);
KN207_1_QL2_20120430_DO_1_dark_uM = KN207_1_QL2_20120430_data(:,3);
KN207_1_QL2_20120430_DO_1_dark_temp_degC = KN207_1_QL2_20120430_data(:,4);

% -------------------------------------------------------------------------
% read in data for NA VICE (KN207-3) Process Station #1, 17-19 June 2012
% -------------------------------------------------------------------------

KN207_3_PS1_20120617 = xlsread('/Users/jrcollins/Code/DO_Instruments/PHORCYS/data/processed/KN207_3/KN207_3_Process_Station_1_PHORCYS_2012_06_17_19.xlsx','Working data','A6:I1403');

KN207_3_PS1_20120617_Timestamp_DO = KN207_3_PS1_20120617(:,1);
KN207_3_PS1_20120617_Timestamp_DO = x2mdate(KN207_3_PS1_20120617_Timestamp_DO,1);
KN207_3_PS1_20120617_timefrac = KN207_3_PS1_20120617(:,2);
KN207_3_PS1_20120617_DO_1_dark_uM_raw = KN207_3_PS1_20120617(:,3);
KN207_3_PS1_20120617_DO_1_dark_sat = KN207_3_PS1_20120617(:,5);
KN207_3_PS1_20120617_DO_1_dark_temp_degC = KN207_3_PS1_20120617(:,6);
KN207_3_PS1_20120617_DO_2_light_uM_raw = KN207_3_PS1_20120617(:,7);
KN207_3_PS1_20120617_DO_2_light_sat = KN207_3_PS1_20120617(:,8);
KN207_3_PS1_20120617_DO_2_light_temp_degC = KN207_3_PS1_20120617(:,9);

% adjust DO values for salinity, see Aanderaa 4330/4835 TD269 manual p. 32

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

% -------------------------------------------------------------------------
% read in data for NA VICE (KN207-3) Process Station #2, 23-27 June 2012,
% -------------------------------------------------------------------------

KN207_3_PS2_20120623 = xlsread('/Users/jrcollins/Code/DO_Instruments/PHORCYS/data/processed/KN207_3/PHORCYS_2012_06_23_27_working.xlsx','Working data','A6:K2930');

KN207_3_PS2_20120623_Timestamp_DO = KN207_3_PS2_20120623(:,1);
KN207_3_PS2_20120623_Timestamp_DO = x2mdate(KN207_3_PS2_20120623_Timestamp_DO,1);
KN207_3_PS2_20120623_timefrac = KN207_3_PS2_20120623(:,2);
KN207_3_PS2_20120623_DO_1_dark_uM_raw = KN207_3_PS2_20120623(:,3);
KN207_3_PS2_20120623_DO_1_dark_sat = KN207_3_PS2_20120623(:,5);
KN207_3_PS2_20120623_DO_1_dark_temp_degC = KN207_3_PS2_20120623(:,6);
KN207_3_PS2_20120623_DO_2_light_uM_raw = KN207_3_PS2_20120623(:,7);
KN207_3_PS2_20120623_DO_2_light_sat = KN207_3_PS2_20120623(:,10);
KN207_3_PS2_20120623_DO_2_light_temp_degC = KN207_3_PS2_20120623(:,11);

% adjust DO values for salinity, see Aanderaa 4330/4835 TD269 manual p. 32

KN207_3_PS2_20120623_DO_1_dark_uM = zeros(length(KN207_3_PS2_20120623_DO_1_dark_uM_raw),1);

for i=1:length(KN207_3_PS2_20120623_DO_1_dark_uM_raw)
    closest_sal = find(KN207_3_mettime>=KN207_3_PS2_20120623_Timestamp_DO(i),1,'first');
    sal = KN207_3_sal(closest_sal);
    T_s = log((298.15-KN207_3_PS2_20120623_DO_1_dark_temp_degC)./(273.15+KN207_3_PS2_20120623_DO_1_dark_temp_degC));
    KN207_3_PS2_20120623_DO_1_dark_uM=KN207_3_PS2_20120623_DO_1_dark_uM_raw.*exp(sal*(B_0+B_1*T_s+B_2*(T_s).^2+B_3*(T_s).^3)+C_0*sal.^2);
end

KN207_3_PS2_20120623_DO_2_light_uM = zeros(length(KN207_3_PS2_20120623_DO_2_light_uM_raw),1);

for i=1:length(KN207_3_PS2_20120623_DO_2_light_uM_raw)
    closest_sal = find(KN207_3_mettime>=KN207_3_PS2_20120623_Timestamp_DO(i),1,'first');
    sal = KN207_3_sal(closest_sal);
    T_s = log((298.15-KN207_3_PS2_20120623_DO_2_light_temp_degC)./(273.15+KN207_3_PS2_20120623_DO_2_light_temp_degC));
    KN207_3_PS2_20120623_DO_2_light_uM=KN207_3_PS2_20120623_DO_2_light_uM_raw.*exp(sal*(B_0+B_1*T_s+B_2*(T_s).^2+B_3*(T_s).^3)+C_0*sal.^2);
end

% normalize light bottle data to calibrated dark bottle data

darklightdiff=KN207_3_PS2_20120623_DO_1_dark_uM-KN207_3_PS2_20120623_DO_2_light_uM;
meandarklightdiff=mean(darklightdiff(96:260,:));

% calibrate dark data to Winklers, using the calibration from the 19 Jun
% deployment at PS#1 since we didn't get Winklers from this station

Wink_diff = Wink_diff_19Jun;
KN207_3_PS2_20120623_DO_1_dark_uM_cal = KN207_3_PS2_20120623_DO_1_dark_uM+Wink_diff;
KN207_3_PS2_20120623_DO_2_light_uM_cal=(KN207_3_PS2_20120623_DO_2_light_uM+meandarklightdiff)+Wink_diff;

% truncate data set to include only usable data
% in this case, we have different ranges for the light and dark bottle data

% lightrange=(105:2730);
% darkrange=(105:2440);

lightrange=(66:2730);
darkrange=(66:2825);

KN207_3_PS2_20120623_Timestamp_DO_dark = KN207_3_PS2_20120623_Timestamp_DO(darkrange,:);
KN207_3_PS2_20120623_timefrac_dark = KN207_3_PS2_20120623_timefrac(darkrange,:);
KN207_3_PS2_20120623_Timestamp_DO_light = KN207_3_PS2_20120623_Timestamp_DO(lightrange,:);
KN207_3_PS2_20120623_timefrac_light = KN207_3_PS2_20120623_timefrac(lightrange,:);
KN207_3_PS2_20120623_DO_1_dark_uM = KN207_3_PS2_20120623_DO_1_dark_uM_cal(darkrange,:);
KN207_3_PS2_20120623_DO_1_dark_temp_degC = KN207_3_PS2_20120623_DO_1_dark_temp_degC(darkrange,:);
KN207_3_PS2_20120623_DO_2_light_uM = KN207_3_PS2_20120623_DO_2_light_uM_cal(lightrange,:);
KN207_3_PS2_20120623_DO_2_light_temp_degC = KN207_3_PS2_20120623_DO_2_light_temp_degC(lightrange,:);

% -------------------------------------------------------------------------
% read in data for NA VICE (KN207-3) Process Station #4, 7-11 July 2012
% -------------------------------------------------------------------------

% rad data not loaded -- clear bottle did not close

KN207_3_PS4_20120711 = xlsread('/Users/jrcollins/Code/DO_Instruments/PHORCYS/data/processed/KN207_3/PHORCYS_2012_07_07_11_working.xlsx','Working data','A6:M2943');

KN207_3_PS4_20120711_Timestamp_DO = KN207_3_PS4_20120711(:,5);
KN207_3_PS4_20120711_Timestamp_DO = x2mdate(KN207_3_PS4_20120711_Timestamp_DO,1);
KN207_3_PS4_20120711_timefrac = KN207_3_PS4_20120711(:,6);
KN207_3_PS4_20120711_DO_1_dark_uM_raw = KN207_3_PS4_20120711(:,7);
KN207_3_PS4_20120711_DO_1_dark_sat = KN207_3_PS4_20120711(:,9);
KN207_3_PS4_20120711_DO_1_dark_temp_degC = KN207_3_PS4_20120711(:,10);

% adjust DO values for salinity, see Aanderaa 4330/4835 TD269 manual p. 32

KN207_3_PS4_20120711_DO_1_dark_uM = zeros(length(KN207_3_PS4_20120711_DO_1_dark_uM_raw),1);

for i=1:length(KN207_3_PS4_20120711_DO_1_dark_uM_raw)
    closest_sal = find(KN207_3_mettime>=KN207_3_PS4_20120711_Timestamp_DO(i),1,'first');
    sal = KN207_3_sal(closest_sal);
    T_s = log((298.15-KN207_3_PS4_20120711_DO_1_dark_temp_degC)./(273.15+KN207_3_PS4_20120711_DO_1_dark_temp_degC));
    KN207_3_PS4_20120711_DO_1_dark_uM=KN207_3_PS4_20120711_DO_1_dark_uM_raw.*exp(sal*(B_0+B_1*T_s+B_2*(T_s).^2+B_3*(T_s).^3)+C_0*sal.^2);
end

% load Winkler data

KN207_3_PS4_dark_Wink_time=xlsread('/Users/jrcollins/Code/DO_Instruments/PHORCYS/data/processed/KN207_3/PHORCYS_2012_07_07_11_working.xlsx','Winkler data','A11');
KN207_3_PS4_dark_Wink_time = x2mdate(KN207_3_PS4_dark_Wink_time,1);
KN207_3_PS4_dark_Wink_uM=xlsread('/Users/jrcollins/Code/DO_Instruments/PHORCYS/data/processed/KN207_3/PHORCYS_2012_07_07_11_working.xlsx','Winkler data','J11');

% calibrate dark data to Winklers

Wink_caltime = find(KN207_3_PS4_20120711_Timestamp_DO>=KN207_3_PS4_dark_Wink_time,1,'first');
Wink_diff = KN207_3_PS4_dark_Wink_uM-KN207_3_PS4_20120711_DO_1_dark_uM(Wink_caltime);
KN207_3_PS4_20120711_DO_1_dark_uM_cal = KN207_3_PS4_20120711_DO_1_dark_uM+Wink_diff;

% truncate data set to include only usable data

KN207_3_PS4_20120711_data = [KN207_3_PS4_20120711_Timestamp_DO KN207_3_PS4_20120711_timefrac ...
    KN207_3_PS4_20120711_DO_1_dark_uM_cal KN207_3_PS4_20120711_DO_1_dark_temp_degC];

% KN207_3_PS4_20120711_data = KN207_3_PS4_20120711_data(180:2893,:);

KN207_3_PS4_20120711_data = KN207_3_PS4_20120711_data(52:end,:);

KN207_3_PS4_20120711_Timestamp_DO = KN207_3_PS4_20120711_data(:,1);
KN207_3_PS4_20120711_timefrac = KN207_3_PS4_20120711_data(:,2);
KN207_3_PS4_20120711_DO_1_dark_uM = KN207_3_PS4_20120711_data(:,3);
KN207_3_PS4_20120711_DO_1_dark_temp_degC = KN207_3_PS4_20120711_data(:,4);

% % -------------------------------------------------------------------------
% % read in data for Palmer Station PHORCYS deployment, 31 Dec 13 - 3 Jan 14
% % -------------------------------------------------------------------------
% 
% % rad data not loaded -- clear bottle did not close
% 
% Pal_20131231 = xlsread('PHORCYS_Palmer_20131231_20140103.xlsx','PHORCYS_20131231_20140103.TXT','F5:M2243');
% 
% Pal_20131231_Timestamp_DO = Pal_20131231(:,1);
% Pal_20131231_Timestamp_DO = x2mdate(Pal_20131231_Timestamp_DO,1);
% Pal_20131231_timefrac = Pal_20131231(:,2);
% Pal_20131231_DO_1_dark_uM_raw = Pal_20131231(:,4);
% Pal_20131231_DO_1_dark_sat = Pal_20131231(:,5);
% Pal_20131231_DO_1_dark_temp_degC = Pal_20131231(:,8);
% 
% % adjust DO values for salinity, see Aanderaa 4330/4835 TD269 manual p. 32
% 
% Pal_20131231_DO_1_dark_uM = zeros(length(Pal_20131231_DO_1_dark_uM_raw),1);
% 
% for i=1:length(Pal_20131231_DO_1_dark_uM_raw)
% %    closest_sal = find(Pal_sal_time>=Pal_20131231_Timestamp_DO(i),1,'first');
% %     sal = Pal_sal(closest_sal);
% %    manually set salinity in absence of complete Palmer salinity data set
%     sal = 32;
%     T_s = log((298.15-Pal_20131231_DO_1_dark_temp_degC)./(273.15+Pal_20131231_DO_1_dark_temp_degC));
%     Pal_20131231_DO_1_dark_uM=Pal_20131231_DO_1_dark_uM_raw.*exp(sal*(B_0+B_1*T_s+B_2*(T_s).^2+B_3*(T_s).^3)+C_0*sal.^2);
% end
% 
% % truncate data set to include only usable data
% 
% Pal_20131231_data = [Pal_20131231_Timestamp_DO Pal_20131231_timefrac ...
%     Pal_20131231_DO_1_dark_uM Pal_20131231_DO_1_dark_temp_degC];
% 
% Pal_20131231_data = Pal_20131231_data(67:1040,:);
% 
% Pal_20131231_Timestamp_DO = Pal_20131231_data(:,1);
% Pal_20131231_timefrac = Pal_20131231_data(:,2);
% Pal_20131231_DO_1_dark_uM = Pal_20131231_data(:,3);
% Pal_20131231_DO_1_dark_temp_degC = Pal_20131231_data(:,4);

%% ------------- define some variables for plots colors etc --------------------------


clf;
close all;

% define some colors

col_t1=[0 0.75 0.75]; % cyan for transparent bottle
col_t2=[0.0784 0.1686 0.55]; % darker cyan for transparent bottle
col_t3=[0.7412 0.8980 0.9765]; % very light cyan for transparent bottle
col_o1=[1 0 0]; % red for opaque bottle
col_o3=[1 .9 .9]; % very light red for opaque bottle
col_temp=[.8 .8 .8]; % a grey for temperature data

% define some font sizes

fsize_ax=8; % for axes
fsize_text=8; % for text

% define a spacing increment for number of minutes in one day

spaceinc = 60*24;

scrsz = get(0,'ScreenSize'); % define a screen size variable so I can make figures that look decent

%% preallocate a matrix to hold our rate data
% and create matrix with indices to the rate data/metadata for each deployment

KN207_Model_A_PHORCYS_met_rates = zeros(5, 8);

KN207_Model_A_PHORCYS_met_rates(:,1) = [2071 2071 2073 2073 2073];
KN207_Model_A_PHORCYS_met_rates(:,2) = [1 2 1 2 4];

KN207_Model_A_PHORCYS_index = table;

% populate our index (dark data, dark timestamps, light data, light
% timestamps)

KN207_Model_A_PHORCYS_index.dark_uM = ...
{'KN207_1_QL1_20120424_DO_1_dark_uM','KN207_1_QL2_20120430_DO_1_dark_uM','KN207_3_PS1_20120617_DO_1_dark_uM','KN207_3_PS2_20120623_DO_1_dark_uM','KN207_3_PS4_20120711_DO_1_dark_uM'}';

KN207_Model_A_PHORCYS_index.dark_timestamp = ...
{'KN207_1_QL1_20120424_Timestamp_DO','KN207_1_QL2_20120430_Timestamp_DO','KN207_3_PS1_20120617_Timestamp_DO','KN207_3_PS2_20120623_Timestamp_DO_dark','KN207_3_PS4_20120711_Timestamp_DO'}';

KN207_Model_A_PHORCYS_index.light_uM = ...
{'','','KN207_3_PS1_20120617_DO_2_light_uM','KN207_3_PS2_20120623_DO_2_light_uM',''}';

KN207_Model_A_PHORCYS_index.light_timestamp = ...
{'','','KN207_3_PS1_20120617_Timestamp_DO','KN207_3_PS2_20120623_Timestamp_DO_light',''}';

%% data analysis

% approach suggested by Emery and Thomson (2001), Data Analysis Meth. in
% P.O. (see chapters 3 and 5), with additional input from Scott Doney

for i=1:size(KN207_Model_A_PHORCYS_met_rates,1)
    
    % first, look at and prepare data for analysis
    
    % remove first-order trend, assumed to be linear, and store the
    % trend in a new variable for each bottle
    
    if strcmp(KN207_Model_A_PHORCYS_index.light_uM(i),'') ~= 1
        
        % there is light bottle data for this station
        
        % light bottle data
        
        light_uM = eval(char(KN207_Model_A_PHORCYS_index.light_uM(i)));
        light_detrended = detrend(eval(char(KN207_Model_A_PHORCYS_index.light_uM(i))));
        light_trend = eval(char(KN207_Model_A_PHORCYS_index.light_uM(i)))-light_detrended;
        light_timestamp = eval(char(KN207_Model_A_PHORCYS_index.light_timestamp(i)));
        
    end
    
    % dark bottle data
    
    dark_uM = eval(char(KN207_Model_A_PHORCYS_index.dark_uM(i)));
    dark_detrended = detrend(eval(char(KN207_Model_A_PHORCYS_index.dark_uM(i))));
    dark_trend = eval(char(KN207_Model_A_PHORCYS_index.dark_uM(i)))-dark_detrended;
    dark_timestamp = eval(char(KN207_Model_A_PHORCYS_index.dark_timestamp(i)));
    
    % plot full and then and detrended data against time
    
    figure
    hold on
    if strcmp(KN207_Model_A_PHORCYS_index.light_uM(i),'') ~= 1
        plot(light_timestamp,light_uM,'r-')
        plot(light_timestamp,light_trend,'r--','LineWidth',2)
    end
    plot(dark_timestamp,dark_uM,'b.')
    plot(dark_timestamp,dark_trend,'b--','LineWidth',2)
    hold off
    if strcmp(KN207_Model_A_PHORCYS_index.light_uM(i),'') ~= 1
        legend('Light bottle','Light bottle - linear trend','Dark bottle','Dark bottle - linear trend');
    else
        legend('Dark bottle','Dark bottle - linear trend');
    end
    ylabel('Dissolved oxygen (uM)');
    datetick('x','dd mmm yy HH:MM');
    set(gca,'XLim',[dark_timestamp(1) dark_timestamp(length(dark_timestamp))]);
    
    % perform calculations
    
    % ------------------------------------------------------------------------
    % light bottle
    % ------------------------------------------------------------------------
    
    if strcmp(KN207_Model_A_PHORCYS_index.light_uM(i),'') ~= 1
        
        % compute autocorrelations of detrended O2 data out to N/3 using unbiased
        % mode; using xcorr rather than autocorr because I have more control over
        % how the function performs
        
        [Cyy,lags]=xcorr(light_detrended,...
            round(length(light_detrended)/3),'unbiased');
        
        % stem(lags,Cyy) % plot
        
        % compute integral from C(0) out to first zero crossing (approximate
        % time scale of decorrelation)
        
        indc = find(lags==0); % index C(0)
        indzero = find(Cyy(indc:length(lags))<=0,1,'first'); % find zero crossing
        
        Z=trapz(Cyy(indc:indc+indzero)); % compute numerical integral of ACF from
        % C(0) to zero crossing
        
        T_unit=Z/max(Cyy); % calculate integral time scale, E&T eqn. 3.15.16a; this value
        % is in unit increments
        T = T_unit*2; % integral time scale in minutes (data were collected at 1/120 Hz)
        
        % now, calculate effective degrees of freedom, N* (E&T eqn. 3.15.17)
        
        N = length(light_detrended); % no. observations
        dt = 2; % sampling interval in minutes
        
        N_star = (N*dt)/T; % N*, the effective # of degrees of freedom, which should be << N
        
        % next, calculate metabolic rats using simple linear regression of each
        % time series (call linfit.m of Glover et al)
        
        [a,sa,cov,r,del,S] = linfit(light_timestamp...
            ,light_uM,0);
        
        met_rate = a(2); % metabolic rate in umol O2/L/d (from slope)
        
        % compute yhats and SSE
        
        Yhats = a(2)*light_timestamp+a(1);
        SSE = sum((light_uM-Yhats).^2);
        
        % compute a realistic 'adjusted' sigma using SSE and N* rather than N
        
        sig = SSE/(N_star-2); % also have to subtract 2 to account for the df's lost in the actual calculation
        
        % use this 'adjusted' sigma to obtain a realistic estimate of uncertainty to
        % accompany the slope
        
        met_rate_uncert_adj = sqrt(sig*S/del);
        
        % store result appropriately
        
        NCP = [met_rate met_rate_uncert_adj T];
        
        KN207_Model_A_PHORCYS_met_rates(i,3:5) = NCP;
        
    end
    
    % ------------------------------------------------------------------------
    % dark bottle
    % ------------------------------------------------------------------------
    
    % compute autocorrelations of detrended O2 data out to N/3 using unbiased
    % mode; using xcorr rather than autocorr because I have more control over
    % how the function performs
    
    [Cyy,lags]=xcorr(dark_detrended,...
        round(length(dark_detrended)/3),'unbiased');
    
    % stem(lags,Cyy) % plot
    
    % compute integral from C(0) out to first zero crossing (approximate
    % time scale of decorrelation)
    
    indc = find(lags==0); % index C(0)
    indzero = find(Cyy(indc:length(lags))<=0,1,'first'); % find zero crossing
    
    Z=trapz(Cyy(indc:indc+indzero)); % compute numerical integral of ACF from
    % C(0) to zero crossing
    
    T_unit=Z/max(Cyy); % calculate integral time scale, E&T eqn. 3.15.16a; this value
    % is in unit increments
    T = T_unit*2; % integral time scale in minutes (data were collected at 1/120 Hz)
    
    % now, calculate effective degrees of freedom, N* (E&T eqn. 3.15.17)
    
    N = length(dark_detrended); % no. observations
    dt = 2; % sampling interval in minutes
    
    N_star = (N*dt)/T; % N*, the effective # of degrees of freedom, which should be << N
    
    % next, calculate metabolic rats using simple linear regression of each
    % time series (call linfit.m of Glover et al)
    
    [a,sa,cov,r,del,S] = linfit(dark_timestamp...
        ,dark_uM,0);
    
    met_rate = a(2); % metabolic rate in umol O2/L/d (from slope)
    
    % compute yhats and SSE
    
    Yhats = a(2)*dark_timestamp+a(1);
    SSE = sum((dark_uM-Yhats).^2);
    
    % compute a realistic 'adjusted' sigma using SSE and N* rather than N
    
    sig = SSE/(N_star-2); % also have to subtract 2 to account for the df's lost in the actual calculation
    
    % use this 'adjusted' sigma to obtain a realistic estimate of uncertainty to
    % accompany the slope
    
    met_rate_uncert_adj = sqrt(sig*S/del);
    
    % store result appropriately
    
    GR = [-met_rate met_rate_uncert_adj T];
    
    KN207_Model_A_PHORCYS_met_rates(i,6:8) = GR;
    
end