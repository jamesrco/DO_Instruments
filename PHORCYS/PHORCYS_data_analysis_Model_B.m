%% PHORCYS_data_analysis_Model_B.m

% Owner/author: James Collins, MIT-WHOI Joint Program, Woods Hole
% Oceanographic Institution, james.r.collins@aya.yale.edu

% Process data from the PHORCYS Model "B" (magnetic drive repeatedly opens
% and closes chambers) and calculate community metabolic rates

% The Model "B" instrument uses Aanderaa Aqua 4531 optodes; salinity was
% set to default of 0s

% Created 11/11/16 by JRC under MATLAB R2015a; maintained on GitHub at
% https://github.com/jamesrco/DO_Instruments under a GNU General Public
% License v. 3

% All time values are in local time

%% -------------------------------------------------------------------------

% clear & close to prep workspace

clear all;

close all;

% set constants for salinity compensation, see Aanderaa 4531 (TD296) manual
% p. 21

B_0 = -6.24097e-3;
B_1 = -6.93498e-3;
B_2 = -6.90358e-3;
B_3 = -4.29155e-3;
C_0 = -3.11680e-7;

% -------------------------------------------------------------------------
% set up for processing of data from November 2016 Iselin dock deployment
% -------------------------------------------------------------------------

% the chambers were programmed to cycle at 0600 and 1700 local time each
% day; these timepoints roughly corresponded to sunrise and sunset

% data were collected at 1-minute intervals

% a new .4KS data file was written for each time segment; this script reads
% in .csv files that have been created from the raw .4KS files after each
% data string consisting of 4 lines was concatenated into a single line

% note that the transparent (clear) chamber did not close for the first 4
% cycles; the DO measured in the clear chamber during these segments
% mirrors the ambient signal and the NCP calculated for these segments is
% thus meaningless

% preallocate a matrix to hold our rate data and create matrix with indices
% to the data files for each time segment; have 8 segments total

Iselin_Nov16_Model_B_PHORCYS_met_rates = zeros(8, 8);

% data stored in Iselin_Nov16_Model_B_PHORCYS_met_rates:
% segment_start_timestamp_local; segment_end_timestamp_local; NCP (rate,
% uncertainty, T; all in umol O2 per L per day); GR (rate, uncertainty, T)

% note that the segment start time is not the same as the file start time
% since the chamber remains open for some fixed period of time during which
% water is exchanged (can tell this from the binary chamber status
% indicator field)

% populate our index (file names of data files for each time
% segment)

Iselin_Nov16_Model_B_PHORCYS_file_index = ...
{'11071700.csv','11080600.csv','11081700.csv','11090600.csv','11091700.csv',...
'11100600.csv','11101700.csv','11110600.csv'}';

%% cycle through each segment, perform calculations & display plot

% these data do not need to be pegged to the initial Winklers since the
% optodes were freshly factory calibrated

% approach for uncertainty estimation suggested by Emery and Thomson
% (2001), Data Analysis Meth. in P.O. (see chapters 3 and 5), with
% additional input from Scott Doney

for i=1:size(Iselin_Nov16_Model_B_PHORCYS_met_rates,1)

    % get file path, open file
    
    fn = char(strcat('/Users/jrcollins/Code/DO_Instruments/PHORCYS/data/processed/Iselin_WHOI_2016_11/',cellstr(Iselin_Nov16_Model_B_PHORCYS_file_index{i})));
    
    % read in data for this segment
    
    PHORdata = readtable(fn,'Delimiter',',','ReadVariableNames',0)
    
    % parse necessary data into different variables
    
    bottle_status = PHORdata(:,6);
    date = PHORdata(:,2);
    time_local = PHORdata(:,3);
    light_uM_raw = PHORdata(:,10);
    light_T_deg_C = PHORdata(:,12);
    dark_uM_raw = PHORdata(:,18);
    dark_T_deg_C = PHORdata(:,18);
        
    % create a timestamp
    
    timestamp_local = datetime(horzcat(char(table2cell(date)),char(table2cell(time_local))),...
        'InputFormat','MM/d/yyHH:mm:ss')
    
    % adjust 
    
    
    