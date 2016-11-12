%% PHORCYS_data_analysis_Model_B.m

% Owner/author: James Collins, MIT-WHOI Joint Program, Woods Hole
% Oceanographic Institution, james.r.collins@aya.yale.edu

% Process data from the PHORCYS Model "B" (magnetic drive repeatedly opens
% and closes chambers) and calculate community metabolic rates

% The Model "B" instrument uses Aanderaa Aqua 4531 optodes; assumes
% salinity was set to default of 0 when data were collected & that user has
% an external source of good salinity data

% Created 11/11/16 by JRC under MATLAB R2015a; maintained on GitHub at
% https://github.com/jamesrco/DO_Instruments under a GNU General Public
% License v. 3

% Dependencies/required files:

% The function linfit.m, available from
% http://www.whoi.edu/12.747/mfiles.html, or from the GitHub repo
% https://github.com/jamesrco/dependencies-useful-scripts

%% clear & close to prep workspace; set some constants

clear all;

close all;

% set constants for salinity compensation, see Aanderaa 4531 (TD296) manual
% p. 21

B_0 = -6.24097e-3;
B_1 = -6.93498e-3;
B_2 = -6.90358e-3;
B_3 = -4.29155e-3;
C_0 = -3.11680e-7;

%% define some variables for plot colors etc

clf;
close all;

% define some font sizes

fsize_ax=8; % for axes
fsize_text=8; % for text

% define a spacing increment for number of minutes in one day

spaceinc = 60*24;

scrsz = get(0,'ScreenSize'); % define a screen size variable so I can make figures that look decent

%% -------------------------------------------------------------------------
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

Iselin_Nov16_Model_B_PHORCYS_met_rates = zeros(8,6);

% data stored in Iselin_Nov16_Model_B_PHORCYS_met_rates:
% NCP (rate, uncertainty, T; all in umol O2 per L per day);
% GR (rate, uncertainty, T)

% note that the segment start time is not the same as the file start time
% since the chamber remains open for some fixed period of time during which
% water is exchanged (can tell this from the binary chamber status
% indicator field)

% populate our index (file names of data files for each time
% segment)

Iselin_Nov16_Model_B_PHORCYS_file_index = ...
{'11071700.csv','11080600.csv','11081700.csv','11090600.csv','11091700.csv',...
'11100600.csv','11101700.csv','11110600.csv'}';

% set salinity; figure is a mean for the deployment period from nearby
% NEERS station at mouth of Waquoit Bay (WAQM3); Woods Hole trend similar
% based on historical data

sal = 32.1;

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
    
    PHORdata_raw = readtable(fn,'Delimiter',',','ReadVariableNames',0);
    
    % first, truncate to omit data collected at the beginning of each
    % segment while the chamber was still open
    
    PHORdata = PHORdata_raw(PHORdata_raw.Var6==1,:);
    
    % parse necessary data into different variables
    
    bottle_status = table2array(PHORdata(:,6));
    date = PHORdata(:,2);
    time_local = PHORdata(:,3);
    light_uM_raw = table2array(PHORdata(:,10));
    light_T_deg_C = table2array(PHORdata(:,12));
    dark_uM_raw = table2array(PHORdata(:,16));
    dark_T_deg_C = table2array(PHORdata(:,18));
    ambient_uM_raw = table2array(PHORdata(:,22));
    ambient_T_deg_C = table2array(PHORdata(:,24));
    
    % create a timestamp
    % all time values for this deployment were recorded in local time
    
    timestamp_local = datetime(horzcat(char(table2cell(date)),char(table2cell(time_local))),...
        'InputFormat','MM/d/yyHH:mm:ss');
    
    % compensate for salinity
    % formula, see Aanderaa 4531 (TD296) manual, p. 21
    
    T_s_light = log((298.15-light_T_deg_C)./(273.15+light_T_deg_C));
    light_uM=light_uM_raw.*exp(sal*(B_0+B_1*T_s_light+B_2*(T_s_light).^2+B_3*(T_s_light).^3)+C_0*sal.^2);
    
    T_s_dark = log((298.15-dark_T_deg_C)./(273.15+dark_T_deg_C));
    dark_uM=dark_uM_raw.*exp(sal*(B_0+B_1*T_s_dark+B_2*(T_s_dark).^2+B_3*(T_s_dark).^3)+C_0*sal.^2);
    
    T_s_ambient = log((298.15-ambient_T_deg_C)./(273.15+ambient_T_deg_C));
    ambient_uM=ambient_uM_raw.*exp(sal*(B_0+B_1*T_s_ambient+B_2*(T_s_ambient).^2+B_3*(T_s_ambient).^3)+C_0*sal.^2);
    
    % remove first-order trend, assumed to be linear, and store the
    % trend in a new variable for each bottle
    
    light_detrended = detrend(light_uM);
    light_trend = light_uM-light_detrended;
    
    dark_detrended = detrend(dark_uM);
    dark_trend = dark_uM-dark_detrended;

    ambient_detrended = detrend(ambient_uM);
    ambient_trend = ambient_uM-ambient_detrended;
    
    % plot full and then and detrended data against time
    
    figure
    hold on
    plot(timestamp_local,light_uM,'r-')
    plot(timestamp_local,light_trend,'r--','LineWidth',2)
    plot(timestamp_local,dark_uM,'b.')
    plot(timestamp_local,dark_trend,'b--','LineWidth',2)
    plot(timestamp_local,ambient_uM,'g.')
    plot(timestamp_local,ambient_trend,'g--','LineWidth',2)

    hold off
    legend('Light bottle','Light bottle - linear trend',...
        'Dark bottle','Dark bottle - linear trend',...
        'Ambient','Ambient - linear trend');
    ylabel('Dissolved oxygen (uM)');
    datetick('x','dd mmm yy HH:MM');
    set(gca,'XLim',[datenum(timestamp_local(1)) datenum(timestamp_local(length(timestamp_local)))]);
    
    % now, can calculate rates by linear regression & make error estimates
    % using the Emery and Thompson method
    
    % ------------------------------------------------------------------------
    % light bottle
    % ------------------------------------------------------------------------
    
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
    
    [a,sa,cov,r,del,S] = linfit(datenum(timestamp_local)...
        ,light_uM,0);
    
    met_rate = a(2); % metabolic rate in umol O2/L/d (from slope)
    
    % compute yhats and SSE
    
    Yhats = a(2)*datenum(timestamp_local)+a(1);
    SSE = sum((light_uM-Yhats).^2);
    
    % compute a realistic 'adjusted' sigma using SSE and N* rather than N
    
    sig = SSE/(N_star-2); % also have to subtract 2 to account for the df's lost in the actual calculation
    
    % use this 'adjusted' sigma to obtain a realistic estimate of uncertainty to
    % accompany the slope
    
    met_rate_uncert_adj = sqrt(sig*S/del);
    
    % store result appropriately
    
    NCP = [met_rate met_rate_uncert_adj T];
    
    Iselin_Nov16_Model_B_PHORCYS_met_rates(i,1:3) = NCP;
    
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
    
    [a,sa,cov,r,del,S] = linfit(datenum(timestamp_local)...
        ,dark_uM,0);
    
    met_rate = a(2); % metabolic rate in umol O2/L/d (from slope)
    
    % compute yhats and SSE
    
    Yhats = a(2)*datenum(timestamp_local)+a(1);
    SSE = sum((dark_uM-Yhats).^2);
    
    % compute a realistic 'adjusted' sigma using SSE and N* rather than N
    
    sig = SSE/(N_star-2); % also have to subtract 2 to account for the df's lost in the actual calculation
    
    % use this 'adjusted' sigma to obtain a realistic estimate of uncertainty to
    % accompany the slope
    
    met_rate_uncert_adj = sqrt(sig*S/del);
    
    % store result appropriately
    
    GR = [-met_rate met_rate_uncert_adj T];
    
    Iselin_Nov16_Model_B_PHORCYS_met_rates(i,4:6) = GR;
    
end
    