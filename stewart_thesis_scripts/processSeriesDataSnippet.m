% processSeriesData.m
%
% Called by ProcessTagData_csv.m
%
% This script formats Series data correctly and makes figures
% User is asked to select a directory of Series csv files (created by wildlife computers
% Argos Message Decoder from a .prv file), or places the address of the dir
% after function name. Data from files "*-Series.csv" and "*-Argos.csv"
% are collected. Figures are created and saved in the same folder as the excel file (optional).
%
% Outside Functions Called:
%    sunrise_.m; sunset_.m -> these two functions call others (see comments)
%    Made by W. Broenkow given by S. Flora at Moss Landing Marine Laboratories
%    DayNight_boxes.m (A. Booth)
%    importArchFile;% (A. Booth Dec. 2008)
%    tagDailyHistos_Shallow.m (J. Stewart Feb 2011)
%    Vertvel_Tags % by J. Stewart March 2011 
% plotOxygenProxy % by J. Stewart March 2011
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 27-Mar-2010 14:07:02, edited 01-Mar-2011: analyses in GMT;
% PST only for labeling, huge redundant overhaul 16-Mar-2011
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : processSeriesData.m

Mfilename = mfilename;

%% FIRST: define a day based on sunrise and sunset 

DateTagV = datevec(DateTag);

% Get calculated sunrise/sunset times:
if exist('sunrise_.m', 'file') %checks to see if W. Broenkow's files are present
    yr = floor(mean(DateTagV(:,1)));
    mon = floor(mean(DateTagV(:,2)));
    da = floor(mean(DateTagV(:,3))); %not exactly logical but sufficient for to find sunrise/sunset
    
    [sr,azimuth] = sunrise_(yr,mon,da,lat,lon*-1,1); %function by W. Broenkow of MLML
    srisea = (datevec(sr/24));
    srise = datenum(0,0,0,srisea(:,4), srisea(:,5),srisea(:,6)); % this way is correct, 17-Mar-2011
    
    [ss,azimuth] = sunset_(yr,mon,da,lat,lon*-1,1); %function by W. Broenkow of MLML
    sseta = (datevec(ss/24));
    sset = datenum(0,0,0,sseta(:,4), sseta(:,5),sseta(:,6));
    
else %approximate
    disp('for accurate sunrise/sunset times download air_sea toolbox from http://woodshole.er.usgs.gov/operations/sea-mat/index.html')
    disp('or get sunrise_.m function by W. Broenkow of Moss Landing Marine Labs')
    disp('approximate sunrise/sunset hours: mar set=1:37GMT rise= 13:39GMT; oct set = 24:50, rise = 13.35; Jul set =2:24, rise = 12:45')
    if mean(Time_at_DepthData.TatDmon)<=6 | mean(Time_at_DepthData.TatDmon)>=10 %Winter
        sset =1;
        srise = 12;
    else % summer
        sset = 2;
        srise = 13;
    end
end

% day/night filter % from Diving_castDO.m by A. Booth. % this way is correct, 17-Mar-2011
DateTagVtimes = datenum(0,0,0,DateTagV(:,4),DateTagV(:,5),DateTagV(:,6));
iNight = (DateTagVtimes<=srise) & (DateTagVtimes>sset);
iDay = (DateTagVtimes>srise) | (DateTagVtimes<=sset);

% don't use DateTag anymore: split into Day and Night and do separately.
% Otherwise the split of the 24-hour period will get wonky. 
timeN = DateTag; timeN(iDay) = NaN; %make day times NaN
timeD = DateTag; timeD(iNight) = NaN; 

DepthN = DepthTag; DepthN(iDay) = NaN;
DepthD = DepthTag; DepthD(iNight) = NaN;

TempN = TempTag; TempN(iDay) = NaN;
TempD = TempTag; TempD(iNight) = NaN;

tagNumD = str2num(tagNumS);

%% setup time

timeNV = datevec(timeN);
timeDV = datevec(timeD);
[timeNyr, timeNmon, timeNda, timeNhr, timeNmn, timeNsec] = datevec(timeN); % break into its elements 
[timeDyr, timeDmon, timeDda, timeDhr, timeDmn, timeDsec] = datevec(timeD); 

% find the unique days and index in order to associate with the correct month.
DeployDaysN = unique(timeNda(~isnan(timeNda)));
DeployDaysD = unique(timeDda(~isnan(timeDda)));

if length(DeployDaysN) == length(DeployDaysD)
    DeployDays = DeployDaysN; % doesn't matter which
elseif length(DeployDaysN) > length(DeployDaysD)
    DeployDays = DeployDaysN;
else
    DeployDays = DeployDaysD;
end
DeployDates = datenum(repmat(yr,length(DeployDays),1), repmat(mon,length(DeployDays),1), DeployDays, 0,0,0); % the year and month won't change from DateTag. But don't use beyond. 
DeployDatesV = datevec(DeployDates);

% %% sunset indices: define a day from sunset to the following sunset %%
b = iNight>0;
bb = diff(b)>0;
bb = [bb; NaN];
bbb = find(bb>0)

% make a vector of sunset-sunset 24-hour period indices
DEPLOYINDEX = repmat(DeployDays(1), bbb(1)-1, 1); % day 1 (the day of deployment)
for i = 1:length(bbb)
    if i < length(bbb)
        deploytemp = repmat(i+DeployDays(1), bbb(i+1)-bbb(i),1);
    else
        deploytemp = repmat(i+DeployDays(1),(length(b)-bbb(i))+1,1);
    end
    DEPLOYINDEX = [DEPLOYINDEX; deploytemp];
end

if tagNumD == 83052 % don't know how to fix this except to ID the exact tag. Problem: one hour recorded on Nov 5 that was logged as Nov 6 w/o sunrise-sunset analysis.
    DEPLOYINDEX = DEPLOYINDEX-1;
end
    
DEPLOYINDEXuni = unique(DEPLOYINDEX);
DEPLOYINDEXuniN = datenum(repmat(yr, length(DEPLOYINDEXuni),1), repmat(mon, length(DEPLOYINDEXuni),1), DEPLOYINDEXuni,0,0,0);

%for display ONLY, change DateTag from GMT to PST (PST = GMT-8h)
DateTagPST = DateTag - datenum(0, 0, 0, 8, 0, 0);

%%

disp('Completed processSeriesData.m')
% ===== EOF [processSeriesData.m] ======
