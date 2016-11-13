% processPTDData.m
%  
% this code was written specifically for tag 83046_09, which was downloaded without time series. 
% 
% called from ProcessTagData_csv.m 
%  
% Outside Functions Called: 
%   importMK10report.m (A. Booth)
%  
% AUTHOR    : J. Stewart jules32@gmail.com 
% DATE      : 27-Mar-2010 14:21:12, revised 27-Feb 2011. Dates crossing months OK. Day/Night fixed 21-Mar-2011  
% Revision  : 1.00  
% DEVELOPED : 7.9.0.529 (R2009b) OSX 
% FILENAME  : processPTDData.m 

% this runs, but clean it up a lot more! created 25-Mar-2010 J. Stewart
Mfilename = mfilename;

%% identify files

disp(' ')
disp('==> Select folder where *-PDT.csv file is contained')
dir1 = uigetdir;
dir1 = [dir1 '/'];
cd(dir1)

% Identify PDT and Histo files
csvfiles = dir(fullfile(dir1,'*.csv')); %file .csv files

PDTfilefun = @(x) char(regexp(x, '.+-PDTs.csv','match')); %files of interest
PDTfilename = cellfun(PDTfilefun,{csvfiles(:).name}, 'UniformOutput',0); %find files of interest within cell array
PDTfiledex = cellfun(@isempty,PDTfilename);
PDTfilename = char(PDTfilename(~PDTfiledex));

Histfilefun = @(x) char(regexp(x, '.+-Histos.csv','match')); %files of interest
Histfilename = cellfun(Histfilefun,{csvfiles(:).name}, 'UniformOutput',0); %find files of interest within cell array
Histfiledex = cellfun(@isempty,Histfilename);
Histfilename = char(Histfilename(~Histfiledex));

%% Import -Histos.csv

fid=fopen(fullfile(dir1,Histfilename));
line = 1;
while line<5
    tline = fgetl(fid);
    if ~isempty(tline) && strcmp(tline(1),'D'), break, end %get header names
    line = line+1;
end
fclose(fid);
Tline = ['''' regexprep(tline,',',''',''') '''']; %reformat
Histos.HeaderNames = eval(['{' Tline '}']); %make into string array

numcol = length(Histos.HeaderNames);

%import column Data
fmat = '%n %n %s %s %s %n %n %n %s';
for s = 1:numcol-9 %make colheaderFormat line
    fmat = strcat(fmat,' %n');
end
fid=fopen(fullfile(dir1,Histfilename));
Histos.RawData = textscan(fid, fmat,'HeaderLines', line, 'Delimiter',',');
fclose(fid);

datecolh = find(cellfun(@(x) strcmp(x,'Date'),Histos.HeaderNames)==1); % find Date column
Histos.RawData{datecolh} = datenum(Histos.RawData{datecolh});

%find lines after pop off (had to modify when WLC changed format such that TAD and TAT are no longer interdispersed)
LocCol = find(cellfun(@(x) strcmp(x,'LocationQuality'),Histos.HeaderNames)==1); %9 find LocationQuality column
locqual = cell2mat(cellfun(@(x) isempty(x),Histos.RawData{1,LocCol},'UniformOutput',0)); % identify when not on the surface (when no locqual)

%remove surface readings
Histos.Data = cellfun(@(x) x(locqual),Histos.RawData,'UniformOutput',0);

%% Import PDT files

fid=fopen(fullfile(dir1,PDTfilename));
line = 1;
while line<5 %get header names
    tline = fgetl(fid);
    if ~isempty(tline) && strcmp(tline(1),'D'), break, end
    line = line+1;
end
fclose(fid);
Tline = ['''' regexprep(tline,',',''';''') '''']; %reformat
PDTs.HeaderNames = eval(['{' Tline '}']); %make into string array

%import column Data
numcol = length(PDTs.HeaderNames);
fmat = '%n %n %s %s %n %n %s';
for s = 1:numcol-7 %make colheaderFormat line adding on %n as long as headers are
    fmat = strcat(fmat,' %n');
end
fid=fopen(fullfile(dir1,PDTfilename));
PDTs.RawData = textscan(fid, fmat,'HeaderLines', line, 'Delimiter',',');
fclose(fid);

datecolp = find(cellfun(@(x) strcmp(x,'Date'),PDTs.HeaderNames) == 1); % find Date column
PDTs.RawData{datecolp} = datenum(PDTs.RawData{datecolp}); %% with the series data tags, the PDT date is very strange. Emailed WLC on May 3 to ask. 
PDTdateV = datevec(PDTs.RawData{datecolp});

%find pop line
LocCol = find(cellfun(@(x) strcmp(x,'LocationQuality'),PDTs.HeaderNames)==1); %find LocationQuality column
locqual = cell2mat(cellfun(@(x) ~isempty(x),PDTs.RawData{1,LocCol},'UniformOutput',0));
locqual2 = find(locqual==1);
if ~isempty(locqual2)
    PDTpopline = locqual2(1)-1;%end of deployment before start transmitting
    
    latcol = find(cellfun(@(x) strcmp(x,'Latitude'),PDTs.HeaderNames)==1); %find Lat column
    latdata = PDTs.RawData{1,latcol};
    lat = latdata(PDTpopline+1); %getting the first lat
    
    loncol = find(cellfun(@(x) strcmp(x,'Longitude'),PDTs.HeaderNames)==1); %find lon column
    londata = PDTs.RawData{1,loncol};
    lon = londata(PDTpopline+1); %getting the first lon
    
else % tag 83051 didn't return the entire PDT record, stole this from the Series record. 
    PDTpopline = length(locqual);
    lon = -118.498;
    lat = 32.092;
end

%get PTT number as a string
pttcol = find(cellfun(@(x) strcmp(x,'Ptt'),PDTs.HeaderNames)==1); %find ptt column
pttdata = PDTs.RawData{1,pttcol};
tagNumS = num2str(pttdata(1));

%% get bins from program report file

PTT = str2double(tagNumS);
year = PDTdateV(1,1);
disp(['Searching for programmer report file PTT#: ' tagNumS ' and year: ' num2str(year)])
try RPTstruct = importMK10report_JS(PTT, year); 
    TempBins = RPTstruct.TempBins';
    DepthBins = RPTstruct.DepthBins';
    prf=1; %indicates the programmer report files was found
catch
    disp(['Error while trying to import ' tagNumS ' programmer report files'])
    err = lasterror;
    disp(err)
    prf=0; %indicates the programmer report files was not found
end

%% reformat Histo "Time-at-Depth Data" and "Time-at-Temperature Data"

HistTypecol = find(cellfun(@(x) strcmp(x,'HistType'),Histos.HeaderNames)==1); %find HistType column
HistTypefilt = Histos.Data{1,HistTypecol}; %pull Data out of HistType column
TatDdex = find(cellfun(@(x) strcmp(x,'TAD'),HistTypefilt)==1); %index numbers for all Time at Depth rows
TatTdex = find(cellfun(@(x) strcmp(x,'TAT'),HistTypefilt)==1); %index numbers for all Time at Temp rows

%separate depth from temp
Histos.TatD = cellfun(@(x) x(TatDdex),Histos.Data,'UniformOutput',0);
Histos.TatT = cellfun(@(x) x(TatTdex),Histos.Data,'UniformOutput',0);

numbinscol = find(cellfun(@(x) strcmp(x,'NumBins'),Histos.HeaderNames)==1); %find numbins column
NumBinsdata = Histos.RawData{1,numbinscol};
maxNumBins = length(DepthBins); %max number of bins 

sumcol = find(cellfun(@(x) strcmp(x,'Sum'),Histos.HeaderNames)==1); %find sum of of all data for that hour column
HistDepthTot = Histos.TatD{1,sumcol}; % total count for Time at Depth
HistTempTot = Histos.TatT{1,sumcol}; 

bincols = cellfun(@(x) char(regexp(x, '^Bin','match')),Histos.HeaderNames,'UniformOutput',0); %find all depth columns
Bincolslog = cell2mat(cellfun(@(x) ~isempty(x),bincols,'UniformOutput',0)); %logical index of depth cols
Bincols = find(Bincolslog==1); % numerical index of depth cols
Bincols = Bincols(1:maxNumBins); % no need to get NaNs from unused bins

%% reformat PDT Data (date, depth, mean temp)

%crop at pop off
PDTs.RawDataCrop = cellfun(@(x) x(1:PDTpopline),PDTs.RawData,'UniformOutput',0); % crop the data once it's popped off

numbinscol = find(cellfun(@(x) strcmp(x,'NumBins'),PDTs.HeaderNames)==1); %find numbins column
NumBinsdata = PDTs.RawDataCrop{1,numbinscol};
maxNumBins = max(NumBinsdata);

depthcols = cellfun(@(x) char(regexp(x, 'Depth.+','match')),PDTs.HeaderNames,'UniformOutput',0); %find all depth columns
Depthcolslog = cell2mat(cellfun(@(x) ~isempty(x),depthcols,'UniformOutput',0)); %logical index of depth cols
Depthcols = find(Depthcolslog==1); % numerical index of depth cols
Depthcols = Depthcols(1:maxNumBins); % no need to get NaN depth Data from unused bins

%pull raw PDT Data
rawdatadates = PDTs.RawDataCrop{1,datecolp};
Seriesdate = rawdatadates; % this is so it will be able to be saved in saveTagDataAsMat.m

rawdatadepth=[];
rawdatamintemp =[];
rawdatamaxtemp =[];
for dc = 1:length(Depthcols)
    rawdatadepth = [rawdatadepth,PDTs.RawDataCrop{1,Depthcols(dc)}];
    rawdatamintemp = [rawdatamintemp,PDTs.RawDataCrop{1,Depthcols(dc)+1}];
    rawdatamaxtemp = [rawdatamaxtemp,PDTs.RawDataCrop{1,Depthcols(dc)+2}];
end

%identify pop offs that stay at depth for a time
pop1=[];
for w = 6:length(rawdatadepth)
    %find where the diff in depth is zero between n and n+1 AND n and n+5 AND where all other bin columns are empty
    if (diff([rawdatadepth(w-5,1),rawdatadepth(w-4,1)]) == 0) && (diff([rawdatadepth(w-5,2),rawdatadepth(w,2)]) == 0) && ~sum(isfinite(rawdatadepth(w-5,:)),2)
        if (nanmean(diff(rawdatadepth(w:end,2)))==0)
            pop1 = w-6;
            break
        end
    end
end
if ~isempty(pop1)
    rawdatadepth = rawdatadepth(1:pop1);
    rawdatamintemp = rawdatamintemp(1:pop1);
    rawdatamaxtemp = rawdatamaxtemp(1:pop1);
    rawdatadates = rawdatadates(1:pop1);
end

%cycle through rows of RawData and vertically concatenate
PDTs.PDTdate = [];
PDTs.PDTdepth = [];
PDTs.TempMean = [];
PDTs.TempMin = [];
PDTs.TempMax = [];
for rd = 1:length(PDTs.RawDataCrop{1,1})
    PDTs.PDTdate = [PDTs.PDTdate;repmat(rawdatadates(rd),8,1)];
    PDTs.PDTdepth = [PDTs.PDTdepth;(rawdatadepth(rd,:))'];
    PDTs.TempMean = [PDTs.TempMean;nanmean([rawdatamintemp(rd,:);rawdatamaxtemp(rd,:)],1)'];
    PDTs.TempMin = [PDTs.TempMin;(rawdatamintemp(rd,:))'];
    PDTs.TempMax = [PDTs.TempMax;(rawdatamaxtemp(rd,:))'];
end

DepthTag = PDTs.PDTdepth;
DateTag = PDTs.PDTdate;
TempTag = PDTs.TempMean;
TempTagMin = PDTs.TempMin;
TempTagMax = PDTs.TempMax;
Seriesdate = PDTs.PDTdate;

%% Setup Dates (day/night): for Histo and PDT   

TADDate = Histos.TatD{1,datecolh}; %output variable
TATDate = Histos.TatT{1,datecolh}; %output variable
TADDateV = datevec(TADDate);
TATDateV = datevec(TATDate);

% Get calculated sunrise/sunset times:
if exist('sunrise_.m', 'file') %checks to see if W. Broenkow's files are present
    yr = floor(mean(TADDateV(:,1)));
    mon = floor(mean(TADDateV(:,2)));
    da = floor(mean(TADDateV(:,3))); %not exactly logical but sufficient for to find sunrise/sunset
    
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
    if mean(TimeAtDepth.TatDmon)<=6 | mean(TimeAtDepth.TatDmon)>=10 %Winter
        sset =1;
        srise = 12;
    else % summer
        sset = 2;
        srise = 13;
    end
end

clear azimuth

%% Apply Day/Night filter to Histos (TimeAtDepth)

TADDateVtimes = datenum(0,0,0,TADDateV(:,4),TADDateV(:,5),TADDateV(:,6)); % TimeAtDepth
TATDateVtimes = datenum(0,0,0,TATDateV(:,4),TATDateV(:,5),TATDateV(:,6)); % TimeAtTemperature

% define filter for day for Depth and Temp
iDayDepth = (TADDateVtimes<=sset) | (TADDateVtimes>=srise);
iDayTemp = (TATDateVtimes<=sset) | (TATDateVtimes>=srise);
iNightDepth = (TADDateVtimes<=srise) & (TADDateVtimes>sset);
iNightTemp = (TATDateVtimes<=srise) & (TATDateVtimes>sset);

% day and night: TIME (Depth)
TADDateNight = TADDate; TADDatenight(iDayDepth) = NaN;
TADDateDay = TADDate; TADDateday(iNightDepth) = NaN;
% day and night: TIME (Temp)
TATDateNight = TATDate; TATDatenight(iDayTemp) = NaN;
TATDateDay = TATDate; TATDateday(iNightTemp) = NaN;
 
% sum all columns and get total percentages
% day and night: DEPTH
HistDepthTotD = sum(HistDepthTot(iDayDepth));  
HistDepthTotN = sum(HistDepthTot(~iDayDepth));
% day and night: TEMP
HistTempTotD = sum(HistTempTot(iDayTemp));
HistTempTotN = sum(HistTempTot(~iDayTemp));

TAD.TatDBins = [];
TAT.TatTBins = [];
for bn = 1:length(Bincols)
    TAD.TatDBins(:,end+1) = Histos.TatD{1,Bincols(bn)}; %output variable
    TAT.TatTBins(:,end+1) = Histos.TatT{1,Bincols(bn)}; %output variable
end

%get percent of time in each bin, day and night
HistDepthBinsD = (nansum(TAD.TatDBins(iDayDepth,:),1)/HistDepthTotD)*100; % DepthDay
HistDepthBinsN = (nansum(TAD.TatDBins(~iDayDepth,:),1)/HistDepthTotN)*100; % DepthNight
HistTempBinsD = (nansum(TAT.TatTBins(iDayTemp,:),1)/HistTempTotD)*100; % TempDay
HistTempBinsN = (nansum(TAT.TatTBins(~iDayTemp,:),1)/HistTempTotN)*100; % TempNight

% make vectors: TIME (Depth)
TADDateNightV = datevec(TADDateNight);
TADDateDayV = datevec(TADDateDay);
[TADNyr, TADNmon, TADNda, TADNhr, TADNmn, TADNsec] = datevec(TADDateNight); % break into its elements 
[TADDyr, TADDmon, TADDda, TADDhr, TADDmn, TADDsec] = datevec(TADDateDay); 

% make vectors: TIME (Temp)
TATDateNightV = datevec(TATDateNight);
TATDateDayV = datevec(TATDateDay);
[TATNyr, TATNmon, TATNda, TATNhr, TATNmn, TATNsec] = datevec(TATDateNight); % break into its elements 
[TATDyr, TATDmon, TATDda, TATDhr, TATDmn, TATDsec] = datevec(TATDateDay);

% find the unique days and index in order to associate with the correct month. Just for one is fine
DeployDaysN = unique(TADNda(~isnan(TADNda)));
DeployDaysD = unique(TADNda(~isnan(TADNda)));

if length(DeployDaysN) == length(DeployDaysD)
    DeployDays = DeployDaysN; % doesn't matter which
elseif length(DeployDaysN) > length(DeployDaysD)
    DeployDays = DeployDaysN;
else
    DeployDays = DeployDaysD;
end
DeployDates = datenum(repmat(TATNyr(1),length(DeployDays),1), repmat(TATNmon(1),length(DeployDays),1), DeployDays, 0,0,0); % the year and month won't change from Seriesdate. But don't use beyond. 

% %% sunset indices: define a day from sunset to the following sunset %%
% DEPTH
b = iDayDepth<1;
bb = diff(b)>0;
bb = [bb; NaN];
bbb = find(bb>0)

% make a vector of sunset-sunset 24-hour period indices
deployINDEXD = repmat(DeployDays(1), bbb(1)-1, 1); % day 1 (the day of deployment)
for i = 1:length(bbb)
    if i < length(bbb)
        deploytemp = repmat(i+DeployDays(1), bbb(i+1)-bbb(i),1);
    else
        deploytemp = repmat(i+DeployDays(1),(length(b)-bbb(i))+1,1);
    end
    deployINDEXD = [deployINDEXD; deploytemp];
end

% TEMP
c = iDayTemp<1;
cc = diff(c)>0;
cc = [cc; NaN];
ccc = find(cc>0)
% make a vector of sunset-sunset 24-hour period indices
deployINDEXT = repmat(DeployDays(1), ccc(1)-1, 1); % day 1 (the day of deployment)
for i = 1:length(ccc)
    if i < length(ccc)
        deploytemp = repmat(i+DeployDays(1), ccc(i+1)-ccc(i),1);
    else
        deploytemp = repmat(i+DeployDays(1),(length(c)-ccc(i))+1,1);
    end
    deployINDEXT = [deployINDEXT; deploytemp];
end
    
deployINDEXTuni = unique(deployINDEXT);
deployINDEXTuniN = datenum(repmat(TATNyr(1), length(deployINDEXTuni),1), repmat(TATNmon(1), length(deployINDEXTuni),1), deployINDEXTuni,0,0,0);

%for display ONLY, change Seriesdate from GMT to PST (PST = GMT-8h)
TADdatePST = TADDate - datenum(0, 0, 0, 8, 0, 0);

%% Apply day/night filter: PDT data  

DateTagV = datevec(DateTag);
DateTagVtimes = datenum(0,0,0,DateTagV(:,4),DateTagV(:,5),DateTagV(:,6));

iNight = (DateTagVtimes<=srise) & (DateTagVtimes>sset);
iDay = (DateTagVtimes>srise) | (DateTagVtimes<=sset);

timeN = DateTag; timeN(iDay) = NaN; % make daytimes Nan 
timeD = DateTag; timeD(iNight) = NaN; 

DepthN = DepthTag; DepthN(iDay) = NaN;
DepthD = DepthTag; DepthD(iNight) = NaN;

TempN = TempTag; TempN(iDay) = NaN;
TempD = TempTag; TempD(iNight) = NaN;

timeNV = datevec(timeN);
timeDV = datevec(timeD);
[timeNyr, timeNmon, timeNda, timeNhr, timeNmn, timeNsec] = datevec(timeN); % break into its elements 
[timeDyr, timeDmon, timeDda, timeDhr, timeDmn, timeDsec] = datevec(timeD); 

% find the unique days and index in order to associate with the correct month.
DeployDaysN = unique(timeNda(~isnan(timeNda)));
DeployDaysD = unique(timeDda(~isnan(timeDda)));

% once figure out vertical migrations, change this (See processSeriesData).
% until then:
buflab = [];

if length(DeployDaysN) == length(DeployDaysD)
    DeployDays = DeployDaysN; % doesn't matter which
elseif length(DeployDaysN) > length(DeployDaysD)
    DeployDays = DeployDaysN;
else
    DeployDays = DeployDaysD;
end

% total time. Not per day. 
DeployTot_da = max(DateTag)-min(DateTag);
DeployTot_hr = DeployTot_da*24;

% %% sunset indices: define a day from sunset to the following sunset %%
b = iNight>0; 
bb = diff(b)>0;
bb = [bb; NaN];
bbb = find(bb>0)

% make a vector of sunset-sunset 24-hour period indices
DEPLOYINDEX = repmat(DeployDays(1), bbb(1)-1, 1); % day 1 (the day of deployment)
DEPLOYINDEXN = datenum(repmat(DateTagV(1,1), bbb(1)-1,1), repmat(DateTagV(1,2), bbb(1)-1,1), DateTagV(1,3),0,0,0);
for i = 1:length(bbb)
    if i < length(bbb)
        deploytemp = repmat(i+DeployDays(1), bbb(i+1)-bbb(i),1);
        DEPLOYINDEXNb = datenum(repmat(yr, bbb(i+1)-bbb(i),1), repmat(mon+1, bbb(i+1)-bbb(i),1), i,0,0,0);
    else
        deploytemp = repmat(i+DeployDays(1),(length(b)-bbb(i))+1,1);
        DEPLOYINDEXNb = datenum(repmat(yr, (length(b)-bbb(i))+1,1), repmat(mon+1, (length(b)-bbb(i))+1,1), i,0,0,0);
    end
    DEPLOYINDEX = [DEPLOYINDEX; deploytemp];
    DEPLOYINDEXN = [DEPLOYINDEXN; DEPLOYINDEXNb];
end
    
% DeployDays replaces DEPLOYINDEXuni
% DeployDates replaces DEPLOYINDEXuniN
% DEPLOYINDEXuni = unique(DEPLOYINDEX);
if str2num(tagNumS) == 83046 % I don't know how else to deal with this since it crosses months
    DEPLOYINDEXunia = datenum(DateTagV(1,1), DateTagV(1,2), DateTagV(1,3),0,0,0);
    DeployDaysb = datenum(repmat(yr, length(DeployDays)-1,1), repmat(mon+1, length(DeployDays)-1,1), DeployDays(1:end-1),0,0,0);
    DeployDates = [DEPLOYINDEXunia; DeployDaysb]; % yes, this overwrites it from above, only for this tag. 
else
    dir1 = [dir1 '_PDTfigs/'];
end

%for display ONLY, change Seriesdate from GMT to PST (PST = GMT-8h)
TADdatePST = TADDate - datenum(0, 0, 0, 8, 0, 0);

%% PDT Setup for Histos DAILY figures below: Calc percent time at depth/temp 

dbar = DepthBins(1:end-1); % there was no data in the 1000+ category %DepthBins(end) = DepthBins(end)+1;
tbar = TempBins(1:end-1);

HistDepthBinsNUtrix = [];
HistDepthBinsDUtrix = [];
HistTempBinsNUtrix = [];
HistTempBinsDUtrix = [];
DeployHours = [];

% depth
for i = 1:length(DeployDays)
    Dlog = deployINDEXD == DeployDays(i);
    
    timeDepthN = TADDateNight(Dlog);
    timeDepthD = TADDateDay(Dlog);
    DepthTagU = HistDepthTot(Dlog);
    
    Depth_binnedNU = hist(DepthN(Dlog),dbar); %count within each bin    
    Depth_binnedDU = hist(DepthD(Dlog),dbar); 
    percDNU = ((Depth_binnedNU/sum(Depth_binnedNU))*100)';
    percDDU = ((Depth_binnedDU/sum(Depth_binnedDU))*100)';
 
    HistDepthBinsNUtrix = [HistDepthBinsNUtrix percDNU];
    HistDepthBinsDUtrix = [HistDepthBinsDUtrix percDDU];
   
    DeployHours = [DeployHours (max(DateTag(Dlog))-min(DateTag(Dlog)))*24];  
end

% temp
for i = 1:length(deployINDEXTuni)
    Tlog = deployINDEXT == deployINDEXTuni(i);
    
    timeTempN = TATDateNight(Tlog);
    timeTempD = TATDateDay(Tlog);
    TempTagU = HistTempTot(Tlog(1:length(HistTempTot)));
     
    Temp_binnedNU = hist(TempTagU(isfinite(timeTempN)),tbar); %count within each bin
    percTNU = ((Temp_binnedNU/sum(Temp_binnedNU))*100)';
    Temp_binnedDU = hist(TempTagU(isfinite(timeTempD(1:length(TempTagU)))),tbar); %count within each bin
    percTDU = ((Temp_binnedDU/sum(Temp_binnedDU))*100)';
    
    HistTempBinsNUtrix = [HistTempBinsNUtrix percTNU];
    HistTempBinsDUtrix = [HistTempBinsDUtrix percTDU];
end

%% Export a time at depth and time at temp summary .csv NEEDED??

% FileOut = [tagNumS 'PDT_DayNightBins.csv'];
% ndheadersC = {'DepthBins','TimeAtDepthNight','TimeAtDepthDay','TempBins','TimeAtTempNight','TimeAtTempDay'};
% 
% FID = fopen(FileOut,'wt');
% % print header cell array
% for H=1:length(ndheadersC)-1
%     fprintf(FID,'%s,',ndheadersC{H});
% end
% fprintf(FID,'%s\n',ndheadersC{H+1}); %print last element with new line character
% %print data:
% fprintf(FID,'%d,%2.3f,%2.3f,%d,%2.3f,%2.3f\n',rot90(fliplr([DepthBins,HistDepthBinsN',HistDepthBinsD',TempBins,HistTempBinsN', HistTempBinsD']))); %fprintf prints each line as each column from top to bottom -- this is dumb
% fclose(FID);
% % type(FileOut) %prints file to screen

%% figure set up

% mini setup
% set up numbering for CompareTags
tagNumD = str2num(tagNumS);
switch tagNumD
    case 83046 % 2009
        tagID = '2';
    case 83052
        tagID = '3';
    case 64004
        tagID = '4';
        %         case 83046 % %2008 tag only: only series data
        %             tagID = '1';
    case 83051
        tagID = '5';
    otherwise
        disp('no matching tag!')
end
tagName = [tagID '-' tagNumS,'-',num2str(yr)];

cmapscale = [2 16]; %CA [5 35]; %Gulf

BinsDN = HistDepthBinsN(1:end-1)';
BinsDD = HistDepthBinsD(1:end-1)';
BinsTN = HistTempBinsN(1:end-1)';
BinsTD = HistTempBinsD(1:end-1)';

summdepthday = summary(DepthTag(isfinite(timeD)));
summtempday = summary(TempTag(isfinite(timeN)));
DepthHistotrix = [repmat(tagNumD,length(BinsDD),1) repmat(DateTagV(1,1),length(BinsDD),1) dbar BinsDN BinsDD];
TempHistotrix = [repmat(tagNumD,length(BinsTD),1) repmat(DateTagV(1,1),length(BinsTD),1) tbar BinsTN BinsTD];

%% plots

% plotHistoData % J. Stewart % problems here, March 21
plotPDTData % J. Stewart

%% plot min and max temp/depth. Not working, March 21 2011

% figure
% subplot(1,2,1)
% h(2) = area(DeployDates, Temptrix(:,2)); % max
% hold on
% h(1) = area(DeployDates, Temptrix(:,1)); % min
% set(h(2),'FaceColor',[0.76 0.87 0.78])
% set(h(1),'FaceColor',[1 1 1])
% set(gca,'Layer','top')
% plot(DeployDates, Temptrix(:,1:2), 'k','LineWidth',2)
% set(gca, 'fontsize', 10, 'fontweight', 'bold')
% ylabel('Temperature [\circC]', 'FontSize', 20, 'FontWeight', 'bold')
% set(gca,'XTick', DeployDates)
% datetick('x','mm-dd','keepticks')
% 
% subplot(1,2,2)
% h(2) = area(DeployDates, Depthtrix(:,2)); % max
% hold on
% h(1) = area(DeployDates, Depthtrix(:,1)); % min
% set(gca, 'YDir', 'reverse')
% set(h(2),'FaceColor',[0.7 0.78 1])
% set(h(1),'FaceColor',[1 1 1])
% set(gca,'Layer','top')
% plot(DeployDates, Depthtrix(:,1:2), 'k','LineWidth',2)
% set(gca, 'fontsize', 10, 'fontweight', 'bold')
% ylabel('Depth [m]', 'FontSize',20, 'FontWeight', 'bold')
% set(gca,'XTick', DeployDates)
% datetick('x','mm-dd','keepticks')
% [ax,h(1)] = suplabel('Date ','x'); 
% [ax,h(2)] = suplabel('','y');
% set(h,'FontSize',20, 'FontWeight', 'bold')
% key = {[tagName ' Daily Min and Max Temp and Depth, PDT Data']; [datestr(min(DeployDates),22) ' - '  datestr(max(DeployDates),22)]};
% [h3] = suptitle(key);
% set(h3,'FontSize',24, 'FontWeight', 'bold')
% FigName = [tagNumS 'MinMaxTD.pdf'];
% annotate_JS(Mfilename, gcf, FigName)
% % orient landscape
% % print('-dpdf', FigName);

%%

%time vs depth max/min
% figure
% set(gcf,'Position',[12   213   842   688])
% plot(PDTdata2.PDTdate, PDTdata2.minDepths*-1, ':', 'Color', [.7 .78 1]); %this is missleading since every hour is split up into 8 bins
% hold on
% plot(PDTdata2.date, PDTdata2.maxDepths*-1, ':', 'Color', [0 0 1]); %this is missleading since every hour is split up into 8 bins
% plot(PDTdata2.date, PDTdata2.minDepths*-1, '.', 'Color', [.7 .78 1]);
% plot(PDTdata2.date, PDTdata2.maxDepths*-1, '.', 'Color', [0 0 1]);
% hold off
% title([tagName ' PDT Data'])
% ylabel('Depth [m]')
% xlabel('Time [GMT]')
% datetick('x',6)
%     eval(['print -dtiff ''' tagNumS 'PDT_timeVdepth_maxmin.tiff'''])%save figure to folder

% %Time series
% figure
% set(gcf,'Position',[12   213   842   688],'Color','w')
% plot(Series.Datetime,Series.Depth*-1,'k')
% %     errorbar(Series.Datetime,Series.Depth*-1,Series.Data{1,10}./2,'k')%comes out very messy
% %     hold on
% %     scatter(Series.Datetime,Series.Depth*-1,20,Series.Temp,'filled');
% %     hold off
% if exist('DayNight_boxes.m','file') %if the function is present the figure will have night boxes
%     [x,y] = DayNight_boxes(DateTag,sr/24,ss/24,get(gca,'Ylim'),0); %will plot boxes
%     figure(5); clf('reset')
%     set(gcf,'Position',[12   213   842   688],'Color','w')
%     fill(x,y,[0.8 0.8 0.8], 'EdgeColor','none')
%     hold on
%     plot(Series.Datetime,Series.Depth*-1,'k')
%     scatter(Series.Datetime,Series.Depth*-1,20,Series.Temp,'filled');
%     hold off
% end
% set(gca,'XLim',[min(Series.Datetime) max(Series.Datetime)])
% title([tagName ' Time Series'])
% ylabel('Depth [m]')
% xlabel('Time [GMT]')
% Hcbar = colorbar('location','SouthOutside', 'LineWidth', 1);
% caxis(cmapscale);%Color axis scaling
% xlabel(Hcbar,'Temperature [\circC]');
% datetick('x',6,'keeplimits')

% eval(['print -dtiff ''' tagNumS 'PDT_TimeSeries_temp.tiff'''])%save figure to folder


%% To save in importTagData_csv

DateRaw = DateTag;
DepthRaw = DepthTag;
TempRaw = TempTag;
TagRaw = repmat(tagNumD, length(DateRaw), 1);

Boxplotrix = [repmat(tagNumD,length(DepthD),1) DEPLOYINDEXN DepthN DepthD TempN TempD]; % DepthD and DepthN are made in ProcessSeriesData

%%
  
disp('Completed processPTDData.m') 
% ===== EOF [processPTDData.m] ======  

