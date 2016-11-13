% compareTagsSST.m

% this is called by manageTagData.m  
% 
% this will do analysis for only CA. If want to change to also GOC, must change directory. 
% 
% Cell titles: 
%         load data from all tags
%         make a filter for shallow depths. See compareTagsSST_EachTag for individual analyses
%         
%  
% Outside Functions Called: 
%       sunrise_.m % MLML Broenkow % this part isn't working currently: trying to break shallow into daytime, nighttime. 
%       sunset_.m % MLML Broenkow
%  
% AUTHOR    : J. Stewart jules32@gmail.com 
% DATE      : 10-Aug-2010 23:35:07; modified for individual tags 22-Feb-2011  
% Revision  : 1.00  
% DEVELOPED : 7.9.0.529 (R2009b) OSX 
% FILENAME  : compareTagsSST.m 

Mfilename = mfilename;

TempBinIncr = .5; % probably smallest you'd want to do is .5, check T error on tags
DepthBinIncr = 25;

% ShallowCutoff defined in manageTagData.m

% should I filter out anytime Depth <5m? Where should I do that? Yikes.
% Probably not, they do come to the surface every so often. And if you
% don't trust that, what else don't you trust? 

%% load data from all tags

dirCompare = '/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/CompareTagDir/_CA_compare'; 
cd(dirCompare);

% import rawTD.mat files 
matfiles = dir(fullfile(dirCompare,'*rawTD.mat')); %find .mat files % CURRENTLY DOES NOT INCLUDE 2008 TAG (moved to different folder)

% load in each file and save data in one big matrix. 
DepthAll = []; % all Depths
TempAll = []; % all Temps
DatesAll = []; % all Dates
DeployDatesAll = []; % deployment dates, but actually just each individual day of tagging
TagNames = []; % individual tag names
DataLength = []; % data length for each tag. 
TagID = [];
for i = 1:length(matfiles)
    infile = matfiles(i).name;
    load(infile)
    DepthAll = [DepthAll; DepthRaw]; 
    TempAll = [TempAll; TempRaw];
    DatesAll = [DatesAll; DateRaw];
    DeployDatesAll = [DeployDatesAll; DeployDates];
    TagNames = [TagNames; TagRaw]; 
    DataLength = [DataLength; length(DepthRaw)];
    TagID = [TagID; repmat(i, length(DepthRaw), 1)]; % add i's so can index the tag properly (ie year)
end

TagNamesUni = unique(TagNames(:,1)); % this wouldn't distinguish between 83046_2008 and 83046_2009
TagIDsUni = unique(TagID);
    
clear infile DeployDates DepthTag TempTag i 

%% make a filter for shallow depths

% for all tags
DepthShallowlog = DepthAll <= ShallowCutoff; % ShallowCutoff set above
DepthShallow = DepthAll(DepthShallowlog,1);
DatesAllShallow = DatesAll(DepthShallowlog,1);
TempShallow = TempAll(DepthShallowlog,1);

length(TempAll) % total amount of data counts for all tags
DataLengthShallow = length(TempShallow) % right now not for individual tags: but need to see that. ***TO DO***

% length(TempAll)     % with all 5 tags         29605
% ShallowCutoff = 30; length(TempShallow)       3909    13%
% ShallowCutoff = 40; length(TempShallow)       6758    23%
% ShallowCutoff = 20; length(TempShallow)       1807    06%
% ShallowCutoff = 10; length(TempShallow)       762     03%

%%
% was trying to run this before committee meeting Nov 23 but not working,
% just go with what I have. This was to try to break it into
% daytime/nighttime. 
% % Sunrise and Sunset data % this if right out of ProcessSeriesData.m
% DatesAllV = datevec(DatesAll);
% % Get calculated sunrise/sunset times:
% if exist('sunrise_.m', 'file') %checks to see if W. Broenkow's files are present
%     yr = floor(mean(DatesAllV(:,1)));    
%     mon = floor(mean(DatesAllV(:,2)));
%     da = floor(mean(DatesAllV(:,3))); %not exactly logical but sufficient for to find sunrise/sunset
%     lat = 36.58;
%     lon = 122.05; % don't do it negatively: see input for sunrise_
% 
%     [sr,azimuth] = sunrise_(yr,mon,da,lat,lon*-1,1); %function by W. Broenkow of MLML
%     srise = (sr/24)-floor(sr/24);
%     sriseN = str2double(datestr(sr/24,'HH'));
%     [ss,azimuth] = sunset_(yr,mon,da,lat,lon*-1,1); %function by W. Broenkow of MLML
%     sset = (ss/24)-floor(ss/24);
%     ssetN = str2double(datestr(ss/24,'HH'));
% else %approximate
%     disp('for accurate sunrise/sunset times download air_sea toolbox from http://woodshole.er.usgs.gov/operations/sea-mat/index.html')
%     disp('or get sunrise_.m function by W. Broenkow of Moss Landing Marine Labs')
%     disp('approximate sunrise/sunset hours: mar set=1:37GMT rise= 13:39GMT; oct set = 24:50, rise = 13.35; Jul set =2:24, rise = 12:45')
% end
% clear azimuth
% 
% % time step 
% DatesAlldiff = diff(DatesAll(1:10));
% increment = DatesAlldiff(1);
% incrementsec = datevec(increment);
% timestep_hrs = incrementsec(5)/60 + incrementsec(6)/3600; %%75seconds
% 
% % day/night filter
% % dayfilterlog = (DatesAllV(:,4)<=sset) | (DatesAllV(:,4)>=srise);
% % DepthDay = DepthAll(dayfilterlog); % used to be called dayfilterD instead of DepthDay
% % DepthNight  = DepthAll(~dayfilterlog); %nightfilterD
% % TempDay  = TempAll(dayfilterlog); %nightfilterT 
% % TempNight = TempAll(~dayfilterlog); %dayfilterT   
% % DatesDay = DatesAll(dayfilterlog);
% % DatesNight = DatesAll(~dayfilterlog);
% 
% %$ maybe this needs to ahppen?
% % make filters for day and night % from Diving_castDO.m by A. Booth
% HrTime = DatesAllV(:,4); % seperate date from time on matlab date format
% iNight = find((HrTime<=sriseN) & (HrTime>ssetN));
% iDay = find((HrTime>sriseN) | (HrTime<=ssetN));
% 
% %calculate duration (sec) of climbs and glides during day and night
% DatesN = DatesAll; DatesN(iDay) = NaN; %make day times NaN
% DatesD = DatesAll; DatesD(iNight) = NaN; %make night times NaN
% 
% DepthN = DepthAll; DepthN(iDay) = NaN;
% DepthD = DepthAll; DepthD(iNight) = NaN;
% 
% TempN = TempAll; TempN(iDay) = NaN;
% TempD = TempAll; TempD(iNight) = NaN;
% %$$
% 
% % filter day for shallow filter
% DepthShallowDaylog = DepthD <= ShallowCutoff;
% DepthShallowDay = DepthD(DepthShallowDaylog);
% DatesShallowDay = DatesD(DepthShallowDaylog);
% TempShallowDay = TempD(DepthShallowDaylog);
% 
% %filter night for shallow filter
% DepthShallowNightlog = DepthN(:,1) <= ShallowCutoff;
% DepthShallowNight = DepthN(DepthShallowNightlog);
% DatesShallowNight = DatesN(DepthShallowNightlog);
% TempShallowNight = TempN(DepthShallowNightlog);
% 
% % some stats
% SumDay = length(DepthD) % total amount of data counts for all tags DAY
% DataLengthShallowDay = length(DepthShallowDay) 
% SumNight = length(DepthN) % total amount of data counts for all tags NIGHT
% DataLengthShallowNight = length(DepthShallowNight) 
% 
% % SumDay                                            16623
% % ShallowCutoff = 30; length(DepthShallowDay)       2924    17.6%
% % ShallowCutoff = 40; length(DepthShallowDay)       4807    29%
% % ShallowCutoff = 20; length(DepthShallowDay)       1321    08%
% % ShallowCutoff = 10; length(DepthShallowDay)       597     03.6%
% 
% % SumNight                                            12982
% % ShallowCutoff = 30; length(DepthShallowNight)       985     7.6%
% % ShallowCutoff = 40; length(DepthShallowNight)       1951    15%
% % ShallowCutoff = 20; length(DepthShallowNight)       486     04%
% % ShallowCutoff = 10; length(DepthShallowNight)       165     01%
% 
% figure; hold on
% [x,y] = DayNight_boxes(DatesAll,srise,sset,get(gca,'Ylim'),0);
% set(gcf,'Position',[12   213   842   688])
% fill(x,y,[0.8 0.8 0.8], 'EdgeColor','none')
% hold on
% plot(DatesAll,DepthAll*-1,'.k')
% % plot(DatesD,DepthD*-1,'.r')
% % plot(DatesN,DepthN*-1, '.')
% hold off
% %         title(tagName,'FontSize',20, 'FontWeight', 'bold')
% %         ylabel('Depth (m)', 'FontSize',16, 'FontWeight', 'bold')
% %         xlabel('Date (PST)', 'FontSize',16, 'FontWeight', 'bold')
% set(gca, 'XLim', [min(DatesAll) max(DatesAll)])
% set(gca, 'fontsize', 12, 'fontweight', 'bold')
% %     set(gca,'XTick', datelabelPSTtrix(2:end)+ datenum(0, 0, 0, 8, 0, 0))
% datetick('x','mm/dd', 'keeplimits','keepticks')
% 
% %visualizing
% figure
% hold on
% plot(DatesAll,DepthAll,'.k')
% plot(DatesD,DepthD,'.r')
% plot(DatesN,DepthN, '.')
% set(gca, 'YDir', 'reverse')
% 
% figure
% plot(DatesShallowNight,DepthShallowNight,'.', 'markersize', 15)

%% Histograms of Shallow. this is modified from processSeriesData.m

% depth bins
binD = floor(min(DepthShallow)):DepthBinIncr:floor(max(DepthShallow)); 
Dvec = nan(length(DepthShallow),1);

for b = 2:length(binD)
    index = find(DepthShallow < binD(b) & DepthShallow > binD(b-1));
    Dvec(index) = binD(b-1);
end

% temp bins
binT = floor(min(TempShallow)):TempBinIncr:floor(max(TempShallow)); 
Tvec = nan(length(DepthShallow),1);

for b = 2:length(binT)
    index = find(TempShallow < binT(b) & TempShallow > binT(b-1));
    Tvec(index) = binT(b-1);
end

clear b index

%% quick stats

ShallowSumm = summary(TempShallow);

%% figures

% quick view of Temp Data
% figure
% hist(TempShallow)%, TempBins) % look at this, and see that I need to split into day/night. 
% figure
% hist(TempShallow, length(binT)) % this just says how many bins

% quick view of Temp data with the specified bins
figure
hist(TempShallow, binT) % this specifies the bins
TempHist = hist(TempShallow, binT);
xlabel('Temperature','FontSize', 14, 'FontWeight', 'bold')
ylabel('Count','FontSize', 14, 'FontWeight', 'bold')
title({['Combined Tag Temperature Histograms: Shallow (<' num2str(ShallowCutoff) 'm)']; ['Median Temp: ' num2str(ShallowSumm(4))]}, 'FontSize', 14, 'FontWeight', 'bold');  
legend(num2str(TagIDsUni))
FigName = ['TempHistosCombinedShallow' num2str(ShallowCutoff) 'm.pdf'];
annotate_JS(Mfilename, gcf, FigName)
% orient landscape % save
% print('-dpdf', FigName);
% % 


figure
hist(TempAll, 1:length(TagIDsUni)) % this specifies the bins
TempHist = hist(TempShallow, binT);
xlabel('Temperature','FontSize', 14, 'FontWeight', 'bold')
ylabel('Data Counts','FontSize', 14, 'FontWeight', 'bold')
title('Count Contributions from Each Individual Tag', 'FontSize', 14, 'FontWeight', 'bold');  
FigName = 'HistCountEachTag.pdf';
annotate_JS(Mfilename, gcf, FigName)
% orient landscape % save
% print('-dpdf', FigName);
% %


%% TimeSeries zoomed in on Shallow area
%%%
%%% go back to previous version of this (on thumbdrive?) that will show
%%% how to make these for each individual tag. Then add this to
% %%% importSeries_csv_JS
% figure
% set(gcf,'Position',[12   213   842   688])
% scatter(SeriesdateShallow,DepthShallow*-1,20,TempShallow,'filled');
% title([tagName ' Seriesdata2'])
% ylabel('Depth [m]')
% xlabel('Time [GMT]')
% Hcbar = colorbar('location','SouthOutside', 'LineWidth', 1);
% %caxis([5 35]);%Color axis scaling % for Mexico
% caxis([4 16]);%Color axis scaling % for CA
% xlabel(Hcbar,'Temperature [\circC]');
% datetick('x',6)
% 
% if exist('DayNight_boxes.m','file') %if the function is present the figure will have night boxes
%     [x,y] = DayNight_boxes(Seriesdate,sr/24,ss/24,get(gca,'Ylim'),0); 
%     figure(200)
%     set(gcf,'Position',[12   213   842   688])
%     fill(x,y,[0.8 0.8 0.8], 'EdgeColor','none')
%     hold on
%     scatter(SeriesdateShallow,DepthShallow*-1,20,TempShallow,'filled');
%     hold off
%     title([tagName ' Time Series Shallow: < ' num2str(ShallowCutoff) 'm'],'FontSize',20, 'FontWeight', 'bold')
%     ylabel('Depth [m]')
%     xlabel('Time [GMT]')
%     datetick('x',6)
%     set(gca, 'XLim', [min(Seriesdate) max(Seriesdate)]);
%     Hcbar = colorbar('location','SouthOutside', 'LineWidth', 1);
%     %caxis([5 35]);%Color axis scaling % for Mexico
%     caxis([4 16]);%Color axis scaling % for CA
%     xlabel(Hcbar,'Temperature [\circC]');
% end
% FigName = [tagNumS 'Series_timeVdepth_tempShallow' num2str(ShallowCutoff) 'm.pdf'];
% annotate_JS(Mfilename, gcf, FigName);
% orient landscape % save
% print('-dpdf', FigName);
% 
% %%
% %%
% %% Histograms of Shallow. this is modified from processSeriesData.m
% 
% %% Make depth and temp bins from Series data
% 
% % depth bins
% DepthMax = max(DepthShallow);
% bins = 10:10:550; % figure out how to automate this based on DepthMax
% Dvec = nan(length(DepthShallow),1);
% 
% for b = 2:length(bins)
%     index = find(DepthShallow < bins(b) & DepthShallow > bins(b-1));
%     Dvec(index) = bins(b-1);
% end
% 
% % temp bins
% TempMin = min(TempShallow);
% TempMax = max(TempShallow);
% bins2 = 5:1:20; % figure out how to automate this based on DepthMax
% Tvec = nan(length(DepthShallow),1);
% 
% for b = 2:length(bins2)
%     index = find(TempShallow < bins2(b) & TempShallow > bins2(b-1));
%     Tvec(index) = bins2(b-1);
% end
% 
% %% Sunrise and Sunset data
%  
% DatesAllV = datevec(SeriesdateShallow);
% 
% % Get calculated sunrise/sunset times:
% if exist('sunrise_.m', 'file') %checks to see if W. Broenkow's files are present
%     yr = floor(mean(DatesAllV(:,1)));    
%     mon = floor(mean(DatesAllV(:,2)));
%     da = floor(mean(DatesAllV(:,3))); %not exactly logical but sufficient for to find sunrise/sunset
% 
%     [sr,azimuth] = sunrise_(yr,mon,da,lat,lon*-1,1); %function by W. Broenkow of MLML
%     srise = str2double(datestr(sr/24,'HH'));
%     
%     [ss,azimuth] = sunset_(yr,mon,da,lat,lon*-1,1); %function by W. Broenkow of MLML
%     sset = str2double(datestr(ss/24,'HH'));
%     
% else %approximate
%     disp('for accurate sunrise/sunset times download air_sea toolbox from http://woodshole.er.usgs.gov/operations/sea-mat/index.html')
%     disp('or get sunrise_.m function by W. Broenkow of Moss Landing Marine Labs')
%     disp('approximate sunrise/sunset hours: mar set=1:37GMT rise= 13:39GMT; oct set = 24:50, rise = 13.35; Jul set =2:24, rise = 12:45')
%     if mean(Time_at_DepthData.TatDmon)<=6 | mean(Time_at_DepthData.TatDmon)>=10 %Winter
%         sset =1;
%         srise = 12;
%     else % summer
%         sset = 2;
%         srise = 13;
%     end
% end
% 
% clear azimuth
% 
% % time step 
% Seriesdatediff = diff(SeriesdateShallow(1:10));
% increment = Seriesdatediff(1);
% incrementsec = datevec(increment);
% timestep_hrs = incrementsec(5)/60 + incrementsec(6)/3600; %%75seconds
% 
% % day/night filter
% dayfilterlog = (DatesAllV(:,4)<=sset) | (DatesAllV(:,4)>=srise);
% dayfilterD = DepthShallow(dayfilterlog);
% nightfilterD = DepthShallow(~dayfilterlog);
% dayfilterT = TempShallow(dayfilterlog);
% nightfilterT = TempShallow(~dayfilterlog);
% 
% % sum all columns and get total percentages 
% SumDay = nansum(dayfilterlog);
% SumNight = nansum(~dayfilterlog);
% 
% % percentages
% binsDay = (SumDay/length(SeriesdateShallow))*100;
% binsNight = (SumNight/length(SeriesdateShallow))*100;
% 
% % make filters for day and night % from Diving_castDO.m by A. Booth
% HrTime = DatesAllV(:,4); % seperate date from time on matlab date format
% iNight = find((HrTime<srise) & (HrTime>sset));
% iDay = find((HrTime>srise) | (HrTime<sset));
% 
% %calculate duration (sec) of climbs and glides during day and night
% timeN = SeriesdateShallow; timeN(iDay) = NaN; %make day times NaN
% timeD = SeriesdateShallow; timeD(iNight) = NaN; %make night times NaN
% 
% DepthN = DepthShallow; DepthN(iDay) = NaN;
% DepthD = DepthShallow; DepthD(iNight) = NaN;
% 
% tagNumD = str2num(tagNumS);
% Boxplotrix = [repmat(tagNumD,length(DepthD),1) SeriesdateShallow DepthN DepthD]; % DepthD and DepthN are made in ProcessSeriesData
% 
% 
% % Calc percent time at depth/temp 
% binsizeD = 5;
% dbar = 0:binsizeD:max(DepthShallow);
% Depth_binnedN = hist(DepthShallow(isfinite(timeN)),dbar); %count within each bin
% percDN = ((Depth_binnedN/sum(Depth_binnedN))*100)';
% Depth_binnedD = hist(DepthShallow(isfinite(timeD)),dbar); %count within each bin
% percDD = ((Depth_binnedD/sum(Depth_binnedD))*100)';
% 
% binsizeT = .5;
% tbar = (floor(min(TempShallow)):binsizeT:max(TempShallow));
% Temp_binnedN = hist(TempShallow(isfinite(timeN)),tbar); %count within each bin
% percTN = ((Temp_binnedN/sum(Temp_binnedN))*100)';
% Temp_binnedD = hist(TempShallow(isfinite(timeD)),tbar); %count within each bin
% percTD = ((Temp_binnedD/sum(Temp_binnedD))*100)';
% 
% %% Setup for DAILY figures below: Calc percent time at depth/temp 
% 
% %for display and analysis, change SeriesdateShallow from GMT to PST (PST = GMT-8h)
% SeriesdatePST = SeriesdateShallow - datenum(0, 0, 0, 8, 0, 0);
% [Seriesyr, Seriesmon, Seriesda, Serieshr, Seriesmn, Seriesec] = datevec(SeriesdatePST);
% SeriesdatePSTV = datevec(SeriesdatePST);
% 
% DeployDays = unique(Seriesda);
% percDNUtrix = [];
% percDDUtrix = [];
% percTNUtrix = [];
% percTDUtrix = [];
% datelabelPSTtrix = [];
% DeployHours = [];
% 
% for i = 1:length(DeployDays)
%     dlog = Seriesda == DeployDays(i);
%     datelabelPST = SeriesdatePST(dlog);
%     timeNU = timeN(dlog);
%     timeDU = timeD(dlog);
%     
%     DepthTagU = DepthShallow(dlog);
%     TempTagU = TempShallow(dlog);
%     
%     Depth_binnedNU = hist(DepthTagU(isfinite(timeNU)),dbar); %count within each bin
%     percDNU = ((Depth_binnedNU/sum(Depth_binnedNU))*100)';
%     Depth_binnedDU = hist(DepthTagU(isfinite(timeDU)),dbar); %count within each bin
%     percDDU = ((Depth_binnedDU/sum(Depth_binnedDU))*100)';
%     
%     Temp_binnedNU = hist(TempTagU(isfinite(timeNU)),tbar); %count within each bin
%     percTNU = ((Temp_binnedNU/sum(Temp_binnedNU))*100)';
%     Temp_binnedDU = hist(TempTagU(isfinite(timeDU)),tbar); %count within each bin
%     percTDU = ((Temp_binnedDU/sum(Temp_binnedDU))*100)';
%     if exist('pttdata', 'var')
%         if sum(isnan(percDDU)) == length(percDDU) %&&  pttdata(1) ~= 83051 % if there are no shallow daytimes.
%             percDDU = percDDU';
%         end
%         if sum(isnan(percTDU)) == length(percTDU) %&&  pttdata(1) ~= 83051% if there are no shallow daytimes.
%             percTDU = percTDU';
%         end
%     end
%     
%     percDNUtrix = [percDNUtrix percDNU];
%     percDDUtrix = [percDDUtrix percDDU];
%     percTNUtrix = [percTNUtrix percTNU];
%     percTDUtrix = [percTDUtrix percTDU];
%     datelabelPSTtrix = [datelabelPSTtrix; datelabelPST(1)];
%     DeployHours = [DeployHours (max(datelabelPST)-min(datelabelPST))*24];
% end
% 
% %% Construct output variables: Seriesdata2,Time_at_TempData,Time_at_DepthData
% 
% Seriesdata2.TagNumber = tagNumS;
% Seriesdata2.month = Seriesmon;
% Seriesdata2.SeriesdatePST = SeriesdatePST;
% Seriesdata2.Seriesdepth = DepthShallow;
% %Seriesdata2.meanTemp = Series.meanTemp; % do I want mean Temp?
% 
% %% Setup for plots
% 
% tagName = [tagNumS,'-',num2str(yr)];
% DeployMonth = unique(Seriesmon);
% DeployDatesV = [repmat(yr,length(DeployDays), 1) repmat(mon,length(DeployDays), 1), DeployDays];
% DeployDates = datenum(DeployDatesV);
% numberofpanels = 6;
% numcols = 3;
% 
% %% Plot Series Histograms 
% 
% summdepthday = summary(DepthShallow(isfinite(timeD)));
% summdepthnight = summary(DepthShallow(isfinite(timeN)));
% summtempday = summary(TempShallow(isfinite(timeD)));
% summtempnight = summary(TempShallow(isfinite(timeN)));
% 
% 
% DepthHistotrix = [repmat(tagNumD,length(percDN),1) repmat(Seriesyr(1),length(percDN),1) dbar' percDN percDD];
% 
% % depth
% figure; 
% H(2) = barh(dbar,percDD,'BarWidth',1);
% hold on
% H(1) = barh(dbar,percDN,'BarWidth',1);
% hold off
% set(H(2),'facecolor',[0.7529    0.7529    0.7529],'EdgeColor',[0.5020    0.5020    0.5020],'LineWidth',1.5)
% set(H(1),'facecolor',[0.5020    0.5020    0.5020],'EdgeColor','k','LineWidth',1.5)
% xlabel('Percent of Time','FontSize', 14, 'FontWeight', 'bold')
% ylabel(['Depth Bin (' num2str(binsizeD) 'm)'],'FontSize', 14, 'FontWeight', 'bold')
% legend(H,'Night','Day','Location','SouthEast')
% set(gca,'YLim',[0,dbar(end)], 'YDir', 'reverse')%[min(dbar), max(dbar)])
% set(gca, 'fontsize', 12, 'fontweight', 'bold')
% box off
% h1 = gca;
% title({[tagName ' Depths, Shallow < ' num2str(ShallowCutoff) 'm']; [datestr(min(SeriesdatePST),22) ' - ' datestr(max(SeriesdatePST),22)];...
%     ['Daytime: med = ' num2str(summdepthday(4), '%3.1f') ', mean = ' num2str(summdepthday(1), '%3.1f') ', n=' num2str(sum(Depth_binnedDU))];...
%     ['Nighttime: med = ' num2str(summdepthnight(4), '%3.1f') ', mean = ' num2str(summdepthnight(1), '%3.1f') ', n=' num2str(sum(Depth_binnedNU))]}, 'fontsize', 14, 'fontweight', 'bold');
% FigName = [tagName '_Series_depthbinsShallow' num2str(ShallowCutoff) 'm.pdf'];
% annotate_JS(Mfilename, gcf, FigName);
% orient landscape % save
% print('-dpdf', FigName);
% 
% % temperature
% figure; 
% H(2) = barh(tbar,percTD,'BarWidth',1);
% hold on
% H(1) = barh(tbar,percTN,'BarWidth',1);
% hold off
% set(H(2),'facecolor',[0.7529    0.7529    0.7529],'EdgeColor',[0.5020    0.5020    0.5020],'LineWidth',1.5)
% set(H(1),'facecolor',[0.5020    0.5020    0.5020],'EdgeColor','k','LineWidth',1.5)
% xlabel('Percent of Time','FontSize', 14, 'FontWeight', 'bold')
% ylabel(['Temperature Bin (' num2str(binsizeT) '\circC)'],'FontSize', 14, 'FontWeight', 'bold')
% legend(H,'Night','Day','Location','SouthEast')
% set(gca,'YLim',[tbar(1),tbar(end)])%[min(dbar), max(dbar)])
% set(gca, 'fontsize', 12, 'fontweight', 'bold')
% box off
% h1 = gca;
% title({[tagName ' Temps, Shallow < ' num2str(ShallowCutoff) 'm']; [datestr(min(SeriesdatePST),22) ' - ' datestr(max(SeriesdatePST),22)];...
%     ['Daytime: med = ' num2str(summtempday(4)) ', mean = ' num2str(summtempday(1), '%3.1f') ', n=' num2str(sum(Temp_binnedDU))];...
%   ['Nighttime: med = ' num2str(summtempnight(4)) ', mean = ' num2str(summtempnight(1), '%3.1f') ', n=' num2str(sum(Temp_binnedNU))]}, 'fontsize', 14, 'fontweight', 'bold');
% FigName = [tagName '_Series_tempbinsShallow' num2str(ShallowCutoff) 'm.pdf'];
% annotate_JS(Mfilename, gcf, FigName);
% orient landscape % save
% print('-dpdf', FigName);
% 
% 
% %% save variables for tag comparison
% 
% % quick fix, change later so that this is more automatic
% tdex = tbar >=11 & tbar <= 15;
% tbar = tbar(tdex);
% percTD = percTD(tdex);
% percTN = percTN(tdex);
% 
% tagNumD = str2num(tagNumS);
% 
% save([CompareTagDir tagName '_HistShallow.mat'] , 'tbar', 'percTD', 'percTN')
% 
% %%
% 
% 
% %%
% 
disp('Completed compareTagsSST.m') 
% ===== EOF [compareTagsSST.m] ======  
