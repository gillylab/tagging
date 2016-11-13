% processSeriesData.m
%
% Called by importTagData_csv.m
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
% plotSeriesData.m % J. Stewart: everything is plotted from here. To clean up this code 24-Mar-2011. 
%    tagDailyHistos_Shallow.m (J. Stewart Feb 2011)
%    Vertvel_Tags % by J. Stewart March 2011 
% plotOxygenProxy % by J. Stewart March 2011
% plotContoursTemp % by J. Stewart March 2011
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 27-Mar-2010 14:07:02, edited 01-Mar-2011: analyses in GMT;
% PST only for labeling, huge redundant overhaul 16-Mar-2011
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : processSeriesData.m

Mfilename = mfilename;

%% mini setup

% set up numbering for CompareTags, and handwrite the timing for the dusk/dawn buffers
tagNumD = str2num(tagNumS);
if TagLocation == 1 % California
    switch tagNumD
        case 64004
            tagID = '4';
        case 83046 % %2008 tag only: only series data
            tagID = '1';
        case 83051
            tagID = '5';
        case 83052
            tagID = '3';
        otherwise
            disp('no matching tag!')
    end
    
elseif TagLocation == 0 % GOC
    switch tagNumD
        case 52869
            tagID = '1';
        case 60970
            tagID = '2';
        case 64006
            tagID = '3';
        case 62007
            tagID = '4';
        case 62009
            tagID = '5';
        case 64004
            tagID = '6';
        case 64006 % this is 2008, the other is 2006. Right now I wont' change this, do that someday. Switch to using the year too: copy code from savetoMat file on _csv. 
            tagID = '7';
        case 83048
            tagID = '8';
        otherwise
            disp('no matching tag!')
    end
end
    
Seriesdata2.TagNumber = tagNumS;
tagName = [tagID '-' tagNumS,'-',datestr(Seriesdate(1), 'yyyy')];
tagNameMat = [tagNumS,'-',datestr(Seriesdate(1), 'yyyy')];

numberofpanels = 6;
numcols = 3;
MaxDepthLim = max(DepthTag)+(0.1*max(DepthTag));

%% FIRST: define a day based on sunrise and sunset 

SeriesdateV = datevec(Seriesdate);

% Get calculated sunrise/sunset times:
if exist('sunrise_.m', 'file') %checks to see if W. Broenkow's files are present
    yr = floor(mean(SeriesdateV(:,1)));
    mon = floor(mean(SeriesdateV(:,2)));
    da = floor(mean(SeriesdateV(:,3))); %not exactly logical but sufficient for to find sunrise/sunset
    
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
SeriesdateVtimes = datenum(0,0,0,SeriesdateV(:,4),SeriesdateV(:,5),SeriesdateV(:,6));

% logical for vertical migrations. Thought about making a separate function dielMigrationID.m after Ashley's findDielMigrations.m, but don't need to be so fancy. 
buffernum1 = 1;
buffernum2 = 2;
buffer1 = datenum(0000,00,00,buffernum1,00,00); %within xx of sunset
buffer2 = datenum(0000,00,00,buffernum2,00,00); %within xx of sunset
iAscent = (SeriesdateVtimes>=0) & (SeriesdateVtimes<(sset+buffer2)); % sunset is right about 01:00, so just do it the hour before to keep it on the same day. 
iDescent = (SeriesdateVtimes>=(srise-buffer1)) & (SeriesdateVtimes<(srise+buffer1));

timeAscent = Seriesdate; timeAscent(~iAscent) = NaN; %make day times NaN
timeDescent = Seriesdate; timeDescent(~iDescent) = NaN; 

DepthAscent = DepthTag; DepthAscent(~iAscent) = NaN;
DepthDescent = DepthTag; DepthDescent(~iDescent) = NaN;

TempAscent = TempTag; TempAscent(~iAscent) = NaN;
TempDescent = TempTag; TempDescent(~iDescent) = NaN;

% don't use Seriesdate anymore: split into Day and Night and do separately. Otherwise the split of the 24-hour period will get wonky.
iNight = (SeriesdateVtimes<=srise) & (SeriesdateVtimes>sset);
iDay = (SeriesdateVtimes>srise) | (SeriesdateVtimes<=sset);
    
if TagDataType == 1
    buflab = [];
    timeN = Seriesdate; timeN(iDay) = NaN; %make day times NaN
    timeD = Seriesdate; timeD(iNight) = NaN;
    
    DepthN = DepthTag; DepthN(iDay) = NaN;
    DepthD = DepthTag; DepthD(iNight) = NaN;
    
    TempN = TempTag; TempN(iDay) = NaN;
    TempD = TempTag; TempD(iNight) = NaN;
    
elseif TagDataType == 0
    dir1 = [dir1 '_noVertMigration/'];
    buflab = ['buff' num2str(buffernum) 'hr'];
    iShallow = iAscent+iDescent+iDay; iShallow = iShallow>0;
    iDeep = iAscent+iDescent+iNight; iDeep = iDeep>0;
    
    timeN = Seriesdate; timeN(iShallow) = NaN; %make alle else NaN
    timeD = Seriesdate; timeD(iDeep) = NaN;
    
    DepthN = DepthTag; DepthN(iShallow) = NaN;
    DepthD = DepthTag; DepthD(iDeep) = NaN;
    
    TempN = TempTag; TempN(iShallow) = NaN;
    TempD = TempTag; TempD(iDeep) = NaN;
    
else error('Input for TagDataType incorrect');
end

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
DeployDates = datenum(repmat(yr,length(DeployDays),1), repmat(mon,length(DeployDays),1), DeployDays, 0,0,0); % the year and month won't change from Seriesdate. But don't use beyond. 

% %% sunset indices: define a day from sunset to the following sunset %%
b = iNight>0;
bb = diff(b)>0;
bb = [bb; NaN];
bbb = find(bb>0)

% make a vector of sunset-sunset 24-hour period indices
DEPLOYINDEX = repmat(DeployDays(1), bbb(1), 1); % day 1 (the day of deployment)
for i = 1:length(bbb)
    if i < length(bbb)
        deploytemp = repmat(i+DeployDays(1), bbb(i+1)-bbb(i),1);
    else
        deploytemp = repmat(i+DeployDays(1),(length(b)-bbb(i)+1),1);
    end
    DEPLOYINDEX = [DEPLOYINDEX; deploytemp];
end

% if tagNumD == 83052 || tagNumD == 83051 % don't know how to fix this except to ID the exact tag. Problem: one hour recorded on Nov 5 that was logged as Nov 6 w/o sunrise-sunset analysis.
    DEPLOYINDEX = DEPLOYINDEX(1:end-1);
% end
    
% DeployDays replaces DEPLOYINDEXuni
% DeployDates replaces DEPLOYINDEXuniN
DEPLOYINDEXN = datenum(repmat(yr,length(DEPLOYINDEX),1), repmat(mon, length(DEPLOYINDEX),1),DEPLOYINDEX,0,0,0); % the time is determined below in DeployHours. don't do it here. April 19 2011. % SeriesdateV(:,4),SeriesdateV(:,5),SeriesdateV(:,6));

%for display ONLY, change Seriesdate from GMT to PST (PST = GMT-8h)
SeriesdatePST = Seriesdate - datenum(0, 0, 0, 8, 0, 0);

%% other setup

Boxplotrix = [repmat(tagNumD,length(DepthD),1) DEPLOYINDEXN DepthN DepthD TempN TempD]; % this will be saved as a .mat file

%% Make histogram bins

% depth bins
binsizeD = 25;
dbar = 0:binsizeD:max(DepthTag);
Depth_binnedN = hist(DepthN,dbar); %count within each bin
Depth_binnedD = hist(DepthD,dbar); 
percDN = ((Depth_binnedN/sum(Depth_binnedN))*100)';
percDD = ((Depth_binnedD/sum(Depth_binnedD))*100)';

% temperature bins
binsizeT = 1;
tbar = (floor(min(TempTag)):0.5:max(TempTag));
Temp_binnedN = hist(TempN,tbar); 
Temp_binnedD = hist(TempD,tbar); 
percTN = ((Temp_binnedN/sum(Temp_binnedN))*100)';
percTD = ((Temp_binnedD/sum(Temp_binnedD))*100)';

%% Setup for DAILY figures below: Calc percent time at depth/temp

% setup for daily matrices of histogram output
DepthNbinned_daily = [];
DepthDbinned_daily = [];
TempNbinned_daily = [];
TempDbinned_daily = [];

DepthNperc_daily = [];
DepthDperc_daily = [];
TempNperc_daily = [];
TempDperc_daily = [];

DeployHours = [];
DepthTagMax = []; % Also get max depth for each day for Hinke model

for i = 1:length(DeployDays)
    dlog = DEPLOYINDEX == DeployDays(i); %i seems to work for GOC, the other way was fine with CA. It's when it spans a month...
    % daily depth bins
    DepthNbinnedi= hist(DepthN(dlog),dbar)'; %count within each bin    
    DepthDbinnedi= hist(DepthD(dlog),dbar)'; 
    DepthNperci = ((DepthNbinnedi/sum(DepthNbinnedi))*100);
    DepthDperci = ((DepthDbinnedi/sum(DepthDbinnedi))*100);
      
    TempNbinnedi = hist(TempN(dlog),tbar)'; 
    TempDbinnedi = hist(TempD(dlog),tbar)'; 
    TempNperci = ((TempNbinnedi/sum(TempNbinnedi))*100);
    TempDperci = ((TempDbinnedi/sum(TempDbinnedi))*100);
    
    if exist('pttdata', 'var')
        if i == 1 && pttdata(1) == 83051 % Don't understand why this is giving a problem
            DepthDperci = DepthDperci';
            TempDperci = TempDperci';
        end
    end
    
    %for troubleshooting
    if size(DepthNperci,2) > size(DepthNperci,1) % if it's a row not column vector
        DepthNperci = DepthNperci';
    end
    if size(DepthDperci,2) > size(DepthDperci,1)
        DepthDperci = DepthDperci';
    end
    if size(TempNperci,2) > size(TempNperci,1)
        TempNperci = TempNperci';
    end
    if size(TempDperci,2) > size(TempDperci,1)
        TempDperci = TempDperci';
    end
    if size(DepthNbinnedi,2) > size(DepthNbinnedi,1) 
        DepthNbinnedi = DepthNbinnedi';
    end
    if size(DepthDbinnedi,2) > size(DepthDbinnedi,1) 
        DepthDbinnedi = DepthDbinnedi';
    end
    if size(TempNbinnedi,2) > size(TempNbinnedi,1) 
        TempNbinnedi = TempNbinnedi';
    end
    if size(TempDbinnedi,2) > size(TempDbinnedi,1) 
        TempDbinnedi = TempDbinnedi';
    end
    
    DepthNbinned_daily = [DepthNbinned_daily DepthNbinnedi]; % ACTUALLY MAYBE NOT NECESSARY
    DepthDbinned_daily = [DepthDbinned_daily DepthDbinnedi];
    TempNbinned_daily = [TempNbinned_daily TempNbinnedi];
    TempDbinned_daily = [TempDbinned_daily TempDbinnedi];
    
    DepthNperc_daily = [DepthNperc_daily DepthNperci];
    DepthDperc_daily = [DepthDperc_daily DepthDperci];
    TempNperc_daily = [TempNperc_daily TempNperci];
    TempDperc_daily = [TempDperc_daily TempDperci];
    
    DeployHours = [DeployHours (max(Seriesdate(dlog))-min(Seriesdate(dlog)))*24]; % this place needs SeriesDate. 
%     i % troubleshoot
%     datestr(min(Seriesdate(dlog)))
%     datestr(max(Seriesdate(dlog)))
    DepthTagMax = [DepthTagMax; max(DepthTag(dlog))];
end

DeployTotLogged_hr = sum(DeployHours);
DeployTotLogged_da = DeployTotLogged_hr/24;
DeployTot_da = max(Seriesdate)-min(Seriesdate);

%% plot

% plot vertical Migrations
plotVertMigration % J. Stewart

% plot all other figs
% plotSeriesData % all regular plotting, put in new script 24-Mar-2011. Also calls plotOxygenProxy % J. Stewart 

% tagDailyHistos_Shallow % J. Stewart % this will also save DepthTagMax. '
% Vertvel_Tags % by J. Stewart. Vertical velocities. % not updated March
% 17.

% plotContoursTemp % J. Stewart, based on A. Carlisle.
% plotContoursTime %

%% To save in importTagData_csv

if TagRecovered
    sampleint = 75;
    lengthM = length(DepthN)/sampleint;
    DateRaw = nan(floor(lengthM),1);
    DepthRaw = nan(floor(lengthM),1);
    TempRaw = nan(floor(lengthM),1);
    BoxplotrixRaw = nan(floor(lengthM),size(Boxplotrix,2));
    
    for i = 1:floor(lengthM)
        DateRaw(i) = Seriesdate(sampleint*i); %%%%%%CHANGE
        DepthRaw(i) = DepthTag(sampleint*i);
        TempRaw(i) = TempTag(sampleint*i);
        BoxplotrixRaw(i,:) = Boxplotrix(sampleint*i,:);
    end
    Boxplotrix = BoxplotrixRaw;
else
    DateRaw = Seriesdate;
    DepthRaw = DepthTag;
    TempRaw = TempTag;
end

TagRaw = repmat(tagNumD, length(DateRaw), 1);
DateRawV = datevec(DateRaw);

clear i sampleint lengthM

%%

disp('Completed processSeriesData.m')
% ===== EOF [processSeriesData.m] ======
