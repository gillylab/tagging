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
iNight = (SeriesdateVtimes<=srise) & (SeriesdateVtimes>sset);
iDay = (SeriesdateVtimes>srise) | (SeriesdateVtimes<=sset);

% don't use Seriesdate anymore: split into Day and Night and do separately.
% Otherwise the split of the 24-hour period will get wonky. 
timeN = Seriesdate; timeN(iDay) = NaN; %make day times NaN
timeD = Seriesdate; timeD(iNight) = NaN; 

DepthN = DepthTag; DepthN(iDay) = NaN;
DepthD = DepthTag; DepthD(iNight) = NaN;

TempN = TempTag; TempN(iDay) = NaN;
TempD = TempTag; TempD(iNight) = NaN;

tagNumD = str2num(tagNumS);
% Boxplotrix = [repmat(tagNumD,length(DepthD),1) $Seriesdate DepthN DepthD TempN TempD]; % this will be saved as a .mat file

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
DeployDatesV = datevec(DeployDates);

a = find(iDay>0);
iDay(a)
aa = find(diff(a) > 1)

a = iDay>0;
aa = diff(a)>0;
aa = [aa; NaN];
aaa = find(aa>0);

DEPLOYINDEX = repmat(DeployDays(1), aaa(1)-1, 1); % day 1 (the day of deployment)
for i = 1:length(aaa)
    if i < length(aaa)
        deploytemp = repmat(i+DeployDays(1), aaa(i+1)-aaa(i),1);
    else
        deploytemp = repmat(i+DeployDays(1),(length(a)-aaa(i))+1,1);
    end
    DEPLOYINDEX = [DEPLOYINDEX; deploytemp];
end
DEPLOYINDEXuni = unique(DEPLOYINDEX);

%for display ONLY, change Seriesdate from GMT to PST (PST = GMT-8h)
SeriesdatePST = Seriesdate - datenum(0, 0, 0, 8, 0, 0);

%% other setup

Seriesdata2.TagNumber = tagNumS;
% Seriesdata2.Seriesdate = Seriesdate; % hasn't this already happende? delete?
% Seriesdata2.Seriesdepth = DepthTag;
%Seriesdata2.meanTemp = Series.meanTemp; % do I want mean Temp?

tagName = [tagNumS,'-',num2str(yr)];
numberofpanels = 6;
numcols = 3;

MaxDepthLim = max(DepthTag)+(0.1*max(DepthTag));

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

for i = 1:length(DeployDates)
    nlog = timeNV(:,3) == DeployDatesV(i,3);
    dlog = timeDV(:,3) == DeployDatesV(i,3);
    
    % daily depth bins
    DepthNbinnedi= hist(DepthN(nlog),dbar)'; %count within each bin    
    DepthDbinnedi= hist(DepthD(dlog),dbar)'; 
    DepthNperci = ((DepthNbinnedi/sum(DepthNbinnedi))*100);
    DepthDperci = ((DepthDbinnedi/sum(DepthDbinnedi))*100);
      
    TempNbinnedi = hist(TempN(nlog),tbar)'; 
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
    
    DepthNbinned_daily = [DepthNbinned_daily DepthNbinnedi]; % ACTUALLY MAYBE NOT NECESSARY
    DepthDbinned_daily = [DepthDbinned_daily DepthDbinnedi];
    TempNbinned_daily = [TempNbinned_daily TempNbinnedi];
    TempDbinned_daily = [TempDbinned_daily TempDbinnedi];
    
    DepthNperc_daily = [DepthNperc_daily DepthNperci];
    DepthDperc_daily = [DepthDperc_daily DepthDperci];
    TempNperc_daily = [TempNperc_daily TempNperci];
    TempDperc_daily = [TempDperc_daily TempDperci];
    
    DeployHours = [DeployHours (max(Seriesdate(dlog))-min(Seriesdate(dlog)))*24]; %$ is this another place I can't let go of Seriesdate?
    DepthTagMax = [DepthTagMax; max(DepthTag(dlog))];
end

%% Plot Series Histograms

summdepthday = summary(DepthTag(isfinite(timeD)));
summdepthnight = summary(DepthTag(isfinite(timeN)));
summtempday = summary(TempTag(isfinite(timeD)));
summtempnight = summary(TempTag(isfinite(timeN)));

DepthHistotrix = [repmat(tagNumD,length(percDN),1) repmat(yr,length(percDN),1) dbar' percDN percDD];
TempHistotrix = [repmat(tagNumD,length(percTN),1) repmat(yr,length(percTN),1) tbar' percTN percTD];

% depth
figure;
H(2) = barh(dbar,percDD,'BarWidth',1);
hold on
H(1) = barh(dbar,percDN,'BarWidth',1);
hold off
set(H(2),'facecolor',[0.7529    0.7529    0.7529],'EdgeColor',[0.5020    0.5020    0.5020],'LineWidth',1.5)
set(H(1),'facecolor',[0.5020    0.5020    0.5020],'EdgeColor','k','LineWidth',1.5)
xlabel('Percent of Time','FontSize', 14, 'FontWeight', 'bold')
ylabel(['Depth Bin (' num2str(binsizeD) 'm)'],'FontSize', 14, 'FontWeight', 'bold')
legend(H,'Night','Day','Location','SouthEast')
set(gca,'YLim',[0,dbar(end)+20], 'YDir', 'reverse')%[min(dbar), max(dbar)])
set(gca, 'fontsize', 12, 'fontweight', 'bold')
box off
h1 = gca;
title({[tagName ' Depths']; ['Day median = ' num2str(summdepthday(4)) ', Night median = ' num2str(summdepthnight(4))];...
    [datestr(min(Seriesdate),22) ' - ' datestr(max(Seriesdate),22) ' (GMT)']}, 'fontsize', 16, 'fontweight', 'bold');
FigName = [tagName '_HistoDepth_bins' num2str(binsizeD) 'm.pdf'];
annotate_JS(Mfilename, gcf, FigName);
orient landscape % save
print('-dpdf', [dir1 FigName]);

% temperature
figure;
H(2) = barh(tbar,percTD,'BarWidth',1);
hold on
H(1) = barh(tbar,percTN,'BarWidth',1);
hold off
set(H(2),'facecolor',[0.7529    0.7529    0.7529],'EdgeColor',[0.5020    0.5020    0.5020],'LineWidth',1.5)
set(H(1),'facecolor',[0.5020    0.5020    0.5020],'EdgeColor','k','LineWidth',1.5)
xlabel('Percent of Time','FontSize', 14, 'FontWeight', 'bold')
ylabel(['Temperature Bin (' num2str(binsizeT) '\circC)'],'FontSize', 14, 'FontWeight', 'bold')
legend(H,'Night','Day','Location','SouthEast')
set(gca,'YLim',[0,tbar(end)+1])%[min(dbar), max(dbar)])
set(gca, 'fontsize', 12, 'fontweight', 'bold')
box off
h1 = gca;
title({[tagName ' Temperatures']; ['Day median = ' num2str(summtempday(4)) ', Night median = ' num2str(summtempnight(4))];...
    ['2.5 percentile = ' num2str(prctile(TempTag(isfinite(timeD)), 2.5)) ', 97.5 percentile = ' num2str(prctile(TempTag(isfinite(timeD)), 97.5))];...
    [datestr(min(Seriesdate),22) ' - ' datestr(max(Seriesdate),22) ' (GMT)']}, 'fontsize', 16, 'fontweight', 'bold');
FigName = [tagName '_HistoTemp_bins' num2str(binsizeT) 'C.pdf'];
annotate_JS(Mfilename, gcf, FigName);
orient landscape % save
print('-dpdf', [dir1 FigName]);

%% Depth histos with OML--should this just go to Tags_CTD?

if OMLhistfig %defined in importTagData_csv
    figure
    H(2) = barh(dbar,percDD,'BarWidth',1);
    hold on
    H(1) = barh(dbar,percDN,'BarWidth',1);
    set(H(2),'facecolor',[0.7529    0.7529    0.7529],'EdgeColor',[0.5020    0.5020    0.5020],'LineWidth',1.5)
    set(H(1),'facecolor',[0.5020    0.5020    0.5020],'EdgeColor','k','LineWidth',1.5)
    xlabel(['Tag % of Time at Depth, ' datestr(min(SeriesdatePST),22) ' - ' datestr(max(SeriesdatePST),22)], ...
        'FontSize', 14, 'FontWeight', 'bold')
    ylabel(['Depth Bin (' num2str(binsizeD) 'm)'],'FontSize', 14, 'FontWeight', 'bold')
    legend(H,'Night','Day','Location','SouthEast')
    h1 = gca;
    set(h1,'YLim',[0,dbar(end)+20], 'YDir', 'reverse')%[min(dbar), max(dbar)])
    set(h1, 'fontsize', 12, 'fontweight', 'bold')
    box off
    title({[tagName ' Depths and MBARI Oxygen  ']; ' '}, 'fontsize', 16, 'fontweight', 'bold');
    plotCTD_Tags % J. Stewart
    FigName = [tagName '_HistoDepth_bins_' num2str(binsizeD) 'm_ROVoxy.pdf'];
    annotate_JS(Mfilename, gcf, FigName);
    dirCTD = '/Users/juliastewart/Thesis/Dg_Tagging/CTDs/';
    %     orient landscape % save
    %     print('-dpdf', [dirCTD FigName]);
end

%% daily histos

% depth
Count = 0;
for i = 1:length(DeployDates)
    Count = Count + 1;
    fignum = floor((Count - 1)./numberofpanels)+1; % set the number for 6-panel figures
    fignum2 = fignum+300;
    panelnum = rem(Count - 1, numberofpanels)+1; % rem takes the remainder after division
    if panelnum == 1
        figure(fignum2); clf;
    end
    subplot(ceil(numberofpanels/numcols), numcols, panelnum); hold on;
    H(2) = barh(dbar,DepthDperc_daily(:,i),'BarWidth',1);
    hold on
    H(1) = barh(dbar,DepthNperc_daily(:,i),'BarWidth',1);
    set(H(2),'facecolor',[0.7529    0.7529    0.7529],'EdgeColor',[0.5020    0.5020    0.5020],'LineWidth',1.5)
    set(H(1),'facecolor',[0.5020    0.5020    0.5020],'EdgeColor','k','LineWidth',1.5)
    title([datestr(DeployDates(i),6) ', ' num2str(ceil(DeployHours(i))) ' hours'], 'FontSize',12, 'FontWeight', 'bold')
    legend(H,'Night','Day','Location','SouthEast')
    set(gca,'YLim',[0, MaxDepthLim], 'YDir', 'reverse')%[min(dbar), max(dbar)])
    set(gca, 'fontsize', 12, 'fontweight', 'bold')
    box off
    h1 = gca;
    if panelnum == 6
        key = {[tagName ' Daily Depth Histograms ' num2str(fignum)]};
        [h3] = suptitle(key);
        [ax,h(1)] = suplabel('Percent of Time ','x'); %('Day of the Year ','x');
        [ax,h(2)] = suplabel(['Depth Bin (' num2str(binsizeD) 'm bins)'],'y');
        set(h3,'FontSize',24, 'FontWeight', 'bold')
        set(h,'FontSize',20, 'FontWeight', 'bold')
        FigName = [tagName '_HistoDepth_bins' num2str(binsizeD) 'm_Daily' num2str(fignum) '.pdf'];
        annotate_JS(Mfilename, gcf, FigName)
        orient landscape % save
        print('-dpdf', [dir1 FigName]);
    end
end
if panelnum < 6 % catch for when there aren't 6 panels on a figure
    key = {[tagName ' Daily Depth Histograms ' num2str(fignum)]};
    [h3] = suptitle(key);
    [ax,h(1)] = suplabel('Percent of Time ','x');
    [ax,h(2)] = suplabel(['Depth Bin (' num2str(binsizeD) 'm bins)'],'y');
    set(h3,'FontSize',24, 'FontWeight', 'bold')
    set(h,'FontSize',20, 'FontWeight', 'bold')
    FigName = [tagName '_HistoDepth_bins' num2str(binsizeD) 'm_Daily' num2str(fignum) '.pdf'];
    annotate_JS(Mfilename, gcf, FigName)
    orient landscape % save
    print('-dpdf', [dir1 FigName]);
end

% temp
Count = 0;
for i = 1:length(DeployDates)
    Count = Count + 1;
    fignum = floor((Count - 1)./numberofpanels)+1; % set the number for 6-panel figures
    fignum2 = fignum+400;
    panelnum = rem(Count - 1, numberofpanels)+1; % rem takes the remainder after division
    if panelnum == 1
        figure(fignum2); clf;
    end
    subplot(ceil(numberofpanels/numcols), numcols, panelnum); hold on;
    H(2) = barh(tbar,TempDperc_daily(:,i),'BarWidth',1);
    hold on
    H(1) = barh(tbar,TempNperc_daily(:,i),'BarWidth',1);
    hold off
    set(H(2),'facecolor',[0.7529    0.7529    0.7529],'EdgeColor',[0.5020    0.5020    0.5020],'LineWidth',1.5)
    set(H(1),'facecolor',[0.5020    0.5020    0.5020],'EdgeColor','k','LineWidth',1.5)
    title([datestr(DeployDates(i),6) ', ' num2str(ceil(DeployHours(i))) ' hours'], 'FontSize',12, 'FontWeight', 'bold')
    legend(H,'Night','Day','Location','SouthEast')
    set(gca,'YLim',[0,tbar(end)+1])
    set(gca, 'fontsize', 12, 'fontweight', 'bold')
    box off
    h1 = gca;
    if panelnum == 6
        key = {[tagName ' Daily Temperature Histograms ' num2str(fignum)]};
        [h3] = suptitle(key);
        [ax,h(1)] = suplabel('Percent of Time ','x'); %('Day of the Year ','x');
        [ax,h(2)] = suplabel('Temperature (\circC bins)','y');
        set(h3,'FontSize',24, 'FontWeight', 'bold')
        set(h,'FontSize',20, 'FontWeight', 'bold')
        FigName = [tagName '_HistoTemp_bins' num2str(binsizeT) 'C_Daily' num2str(fignum) '.pdf'];
        annotate_JS(Mfilename, gcf, FigName)
        orient landscape % save
        print('-dpdf', [dir1 FigName]);
    end
end
if panelnum < 6 % catch for when there aren't 6 panels on a figure
    key = {[tagName ' Daily Temperature Histograms ' num2str(fignum)]};
    [h3] = suptitle(key);
    [ax,h(1)] = suplabel('Percent of Time ','x');
    [ax,h(2)] = suplabel('Temperature (\circC bins)','y');
    set(h3,'FontSize',24, 'FontWeight', 'bold')
    set(h,'FontSize',20, 'FontWeight', 'bold')
    FigName = [tagName '_HistoTemp_bins' num2str(binsizeT) 'C_Daily' num2str(fignum) '.pdf'];
    annotate_JS(Mfilename, gcf, FigName)
    orient landscape % save
    print('-dpdf', [dir1 FigName]);
end

%% daily boxplots

Count = 0;
for i = 1:length(DEPLOYINDEXuni)
    %     nlog = timeNV(:,3) == DeployDatesV(i,3);
    %     dlog = timeDV(:,3) == DeployDatesV(i,3);
    dlog = DEPLOYINDEX == DEPLOYINDEXuni(i);
    DepthNn = DepthN(dlog);
    DepthDd = DepthD(dlog);
    if length(DepthNn) > length(DepthDd) % vectors must be the same length
        DepthDD = [DepthDd; nan((length(DepthNn)-length(DepthDd)),1)];
        DepthNN = DepthNn;
    else
        DepthDD = DepthDd;
        DepthNN = [DepthNn; nan((length(DepthDd)-length(DepthNn)),1)];
    end
    boxplotND = [DepthNN DepthDD];
    Count = Count + 1;
    fignum = floor((Count - 1)./numberofpanels)+1; % set the number for 6-panel figures
    fignum2 = fignum+260;
    panelnum = rem(Count - 1, numberofpanels)+1; % rem takes the remainder after division
    if panelnum == 1
        figure(fignum2); clf;
    end
    subplot(ceil(numberofpanels/numcols), numcols, panelnum); hold on;
    boxplot(boxplotND, 'colors', 'k', 'symbol', '.k', 'outliersize', 5, 'labels', {'Night', 'Day'});
    set(gca, 'YLim', [0 MaxDepthLim],'YDir', 'reverse')
    
    title([datestr(DeployDates(i),6) ', ' num2str(ceil(DeployHours(i))) ' hours'], 'FontSize',12, 'FontWeight', 'bold')
    if panelnum == 6
%         fixfig
        key = {[tagName ' Daily Depth Boxplots ' num2str(fignum)]};
        [h3] = suptitle(key);
        [ax,h(2)] = suplabel('Depth (m)','y');
        set(h3,'FontSize',24, 'FontWeight', 'bold')
        set(h,'FontSize',20, 'FontWeight', 'bold')
        FigName = [tagName '_BoxplotDepthDaily' num2str(fignum) '.pdf'];
        annotate_JS(Mfilename, gcf, FigName)
        orient landscape % save
        print('-dpdf', [dir1 FigName]);
    end
end
if panelnum < 6 % catch for when there aren't 6 panels on a figure
%     fixfig
    key = {[tagName ' Daily Depth Boxplots ' num2str(fignum)]};
    [h3] = suptitle(key);
    [ax,h(2)] = suplabel('Depth (m)','y');
    set(h3,'FontSize',24, 'FontWeight', 'bold')
    set(h,'FontSize',20, 'FontWeight', 'bold')
    FigName = [tagName '_BoxplotDepthDaily' num2str(fignum) '.pdf'];
    annotate_JS(Mfilename, gcf, FigName)
    orient landscape % save
    print('-dpdf', [dir1 FigName]);
end

%%

tagDailyHistos_Shallow % J. Stewart % this will also save DepthTagMax

%% temp vs. depth

% temp vs. depth for whole tag
figure;
set(gcf,'Position',[12   213   842   688])
plot(TempTag, DepthTag, '.', 'Color', [0 0.5020 0]);
ylabel('Depth [m]', 'FontSize',20, 'FontWeight', 'bold')
xlabel('Temperature [\circC]', 'FontSize',20, 'FontWeight', 'bold')
title([tagName ' Temperature v. Depth, SeriesData'], 'FontSize',20, 'FontWeight', 'bold')
set(gca, 'YDir', 'reverse');
FigName = [tagNumS '_Profile_TempDepth.pdf'];
annotate_JS(Mfilename, gcf, FigName)
orient landscape % save
print('-dpdf', FigName);

% temp vs. depth daily subplots
Count = 0;
percx = .6;
percy = .75;

Temptrix = zeros(length(DeployDates), 2);
Depthtrix = zeros(length(DeployDates), 2);

for i = 1:length(DeployDates)
    nlog = timeNV(:,3) == DeployDatesV(i,3);
    dlog = timeDV(:,3) == DeployDatesV(i,3);
    DepthNn = DepthN(nlog);
    DepthDd = DepthD(dlog);
    Count = Count + 1;
    fignum = floor((Count - 1)./numberofpanels)+1; % set the number for 6-panel figures
    fignum2 = fignum+100;
    panelnum = rem(Count - 1, numberofpanels)+1; % rem takes the remainder after division
    if panelnum == 1
        figure(fignum2); clf;
    end
    subplot(ceil(numberofpanels/numcols), numcols, panelnum); hold on;
    plot(TempD(dlog), DepthD(dlog), '.', 'Color', [0 0.5020 0]);
    plot(TempN(nlog), DepthN(nlog), '.', 'Color', [0 0.5020 0]);
    set(gca,'XLim', [min(TempTag) max(TempTag)], 'YLim', [min(DepthTag), max(DepthTag)]);
    text(max(TempTag)*percx, max(DepthTag)*percy, {['Tmin = ' num2str(min(TempTag(dlog))) ' ']; ['Tmax = ' num2str(max(TempTag(dlog))) ' ']; ...
        ['Dmin = ' num2str(min(DepthTag(dlog))) ' ']; ['Dmax = ' num2str(max(DepthTag(dlog))) ' ']}, 'fontsize', 12, 'fontweight', 'bold')
    title([datestr(DeployDates(i),6) ', ' num2str(ceil(DeployHours(i))) ' hours'], 'FontSize',12, 'FontWeight', 'bold') 
    set(gca, 'YDir', 'reverse');
    Temptrix(i,:) = [min(TempTag(dlog)) max(TempTag(dlog))];
    Depthtrix(i,:) = [min(DepthTag(dlog)), max(DepthTag(dlog))];
    if panelnum == 6
        key = {[tagName ' Daily Temperature/Depth Profiles ' num2str(fignum)]};
        [h3] = suptitle(key);
        [ax,h(1)] = suplabel('Depth [m] ','x'); %('Day of the Year ','x');
        [ax,h(2)] = suplabel('Temperature [\circC] ','y');
        set(h3,'FontSize',24, 'FontWeight', 'bold')
        set(h,'FontSize',20, 'FontWeight', 'bold')
        FigName = [tagNumS '_Profile_TempDepth_Daily' num2str(fignum) '.pdf'];
        annotate_JS(Mfilename, gcf, FigName)
        orient landscape % save
        print('-dpdf', FigName);
    end
end
if panelnum < 6 % catch for when there aren't 6 panels on a figure
    key = {[tagName ' Daily Temperature/Depth Profiles ' num2str(fignum)]};
    [h3] = suptitle(key);
    [ax,h(1)] = suplabel('Temperature [\circC] ','x');
    [ax,h(2)] = suplabel('Depth [m] ','y');
    set(h3,'FontSize',24, 'FontWeight', 'bold')
    set(h,'FontSize',20, 'FontWeight', 'bold')
    FigName = [tagNumS '_Profile_TempDepth_Daily' num2str(fignum) '.pdf'];
    annotate_JS(Mfilename, gcf, FigName)
    orient landscape % save
    print('-dpdf', FigName);
end

%% plot min and max temp/depth.

figure
subplot(1,2,1)
h(2) = area(DeployDates, Temptrix(:,2)); % max
hold on
h(1) = area(DeployDates, Temptrix(:,1)); % min
set(h(2),'FaceColor',[0.76 0.87 0.78])
set(h(1),'FaceColor',[1 1 1])
set(gca,'Layer','top')
plot(DeployDates, Temptrix(:,1:2), 'k','LineWidth',2)
set(gca, 'fontsize', 12, 'fontweight', 'bold')
ylabel('Temperature [\circC]', 'FontSize', 20, 'FontWeight', 'bold')
set(gca,'XTick', DeployDates)
datetick('x','mm-dd','keepticks')

subplot(1,2,2)
h(2) = area(DeployDates, Depthtrix(:,2)); % max
hold on
h(1) = area(DeployDates, Depthtrix(:,1)); % min
set(gca, 'YDir', 'reverse')
set(h(2),'FaceColor',[0.7 0.78 1])
set(h(1),'FaceColor',[1 1 1])
set(gca,'Layer','top')
plot(DeployDates, Depthtrix(:,1:2), 'k','LineWidth',2)
set(gca, 'fontsize', 12, 'fontweight', 'bold')
ylabel('Depth [m]', 'FontSize',20, 'FontWeight', 'bold')
set(gca,'XTick', DeployDates)
datetick('x','mm-dd','keepticks')
[ax,h(1)] = suplabel('Date ','x');
[ax,h(2)] = suplabel('','y');
set(h,'FontSize',20, 'FontWeight', 'bold')
key = {[tagName ' Daily Min and Max Temp and Depth, Series Data']; [datestr(min(DeployDates),22) ' - '  datestr(max(DeployDates),22)]};
[h3] = suptitle(key);
set(h3,'FontSize',24, 'FontWeight', 'bold')
FigName = [tagNumS 'MinMaxTD.pdf'];
annotate_JS(Mfilename, gcf, FigName)
orient landscape
print('-dpdf', FigName);

% $fixfig mlc

clear ax

%% datetime vs depth with temp as color
% 
% if ~TagRecovered % because Recovered Tag Time Series takes so long to plot, don't do it automatically!
%     figure(200); clf('reset');
%     set(gcf,'Position',[12   213   842   688])
%     scatter(Seriesdate,DepthTag*-1,20,TempTag,'filled');$
%     title([tagName ' Seriesdata2'])
%     ylabel('Depth (m)')
%     xlabel('Date (PST)')
%     Hcbar = colorbar('location','SouthOutside', 'LineWidth', 1);
%     %caxis([5 35]);%Color axis scaling % for Mexico
%     caxis([4 16]);%Color axis scaling % for CA
%     xlabel(Hcbar,'Temperature (\circC)');
%     datetick('x',6)
%     if exist('DayNight_boxes.m','file') %if the function is present the figure will have night boxes
%         [x,y] = DayNight_boxes(Seriesdate,sr/24,ss/24,get(gca,'Ylim'),0);
%         figure(200); clf('reset');
%         set(gcf,'Position',[12   213   842   688])
%         fill(x,y,[0.8 0.8 0.8], 'EdgeColor','none')
%         hold on
%         scatter(Seriesdate,DepthTag*-1,20,TempTag,'filled');
%         hold off
%         title(tagName,'FontSize',20, 'FontWeight', 'bold')
%         ylabel('Depth (m)', 'FontSize',16, 'FontWeight', 'bold')
%         xlabel('Date (PST)', 'FontSize',16, 'FontWeight', 'bold')
%         set(gca, 'XLim', [min(Seriesdate) max(Seriesdate)])
%         set(gca, 'fontsize', 12, 'fontweight', 'bold')
%         Hcbar = colorbar('location','SouthOutside', 'LineWidth', 1);
%         %caxis([5 35]);%Color axis scaling % for Mexico
%         caxis([4 16]);%Color axis scaling % for CA
%         xlabel(Hcbar,'Temperature (\circC)', 'FontSize',16, 'FontWeight', 'bold');
%         set(gca,'XTick', DeployDates+ datenum(0, 0, 0, 8, 0, 0))
%         datetick('x','mm/dd', 'keeplimits','keepticks')
%     end
%     FigName = [tagNumS 'Series_timeVdepth_temp.pdf'];
%     annotate_JS(Mfilename, gcf, FigName);
%     orient landscape % save
%     print('-dpdf', FigName);
%     
%     plotOxygenProxy % J. Stewart
%     
% end

%% Look at vertical velocity to determine a rate each day.


%% troubleshooting day.night

% figure
% scatter(Seriesdate,DepthTag*-1,20,TempTag,'filled');
% hold on
% plot(Seriesdate(iNight),DepthTag(iNight)*-1, '.k', 'MarkerSize', 15);
% title([tagName ' Seriesdata2'])
% ylabel('Depth (m)')
% xlabel('Date (PST)')
% Hcbar = colorbar('location','SouthOutside', 'LineWidth', 1);
% %caxis([5 35]);%Color axis scaling % for Mexico
% caxis([4 16]);%Color axis scaling % for CA
% xlabel(Hcbar,'Temperature (\circC)');
% datetick('x',6)

    [x,y] = DayNight_boxes(Seriesdate,sr/24,ss/24,get(gca,'Ylim'),0);
%  [x,y] = DayNight_boxes(Seriesdate,datenum(0000,00,00,srise,00,00),datenum(0000,00,00,sset,00,00),get(gca,'Ylim'),0);
    figure;
    set(gcf,'Position',[12   213   842   688])
    fill(x,y,[0.8 0.8 0.8], 'EdgeColor','none')
    hold on
    scatter(Seriesdate(iNight),DepthTag(iNight)*-1,20,TempTag(iNight),'filled');
    scatter(Seriesdate(iDay),DepthTag(iDay)*-1,20,TempTag(iDay),'filled');
%     plot(Seriesdate(iNight),DepthTag(iNight)*-1, '.k', 'MarkerSize', 15);
%     plot(Seriesdate(iDay),DepthTag(iDay)*-1, '.r', 'MarkerSize', 15);
    hold off
    title(tagName,'FontSize',20, 'FontWeight', 'bold')
    ylabel('Depth (m)', 'FontSize',16, 'FontWeight', 'bold')
    xlabel('Date (PST)', 'FontSize',16, 'FontWeight', 'bold')
    set(gca, 'XLim', [min(Seriesdate) max(Seriesdate)])
    set(gca, 'fontsize', 12, 'fontweight', 'bold')
    Hcbar = colorbar('location','SouthOutside', 'LineWidth', 1);
    %caxis([5 35]);%Color axis scaling % for Mexico
    caxis([4 16]);%Color axis scaling % for CA
    xlabel(Hcbar,'Temperature (\circC)', 'FontSize',16, 'FontWeight', 'bold');
    set(gca,'XTick', DeployDates+ datenum(0, 0, 0, 8, 0, 0))
    datetick('x','mm/dd', 'keeplimits','keepticks')
    
%     hold on
%     for i = 1:length(DeployDates)
%         dlog = SeriesdateV(:,3) == DeployDatesV(i,3);
%         plot(Seriesdate(dlog),DepthTag(dlog)*-1, '.b', 'MarkerSize', 15);
%     end

    

%% To save in importTagData_csv

if TagRecovered
    sampleint = 75;
    lengthM = length(DepthN)/sampleint;
    DateRaw = nan(floor(lengthM),1);
    DepthRaw = nan(floor(lengthM),1);
    TempRaw = nan(floor(lengthM),1);
    
    for i = 1:floor(lengthM)
        DateRaw(i) = Seriesdate(sampleint*i); %%%%%%CHANGE
        DepthRaw(i) = DepthTag(sampleint*i);
        TempRaw(i) = TempTag(sampleint*i);
    end
else
    DateRaw = Seriesdate;
    DepthRaw = DepthTag;
    TempRaw = TempTag;
end

TagRaw = repmat(tagNumD, length(DateRaw), 1);
DateRawV = datevec(DateRaw);

clear i sampleint lengthM

%%

Vertvel_Tags % by J. Stewart

%%

disp('Completed processSeriesData.m')
% ===== EOF [processSeriesData.m] ======
