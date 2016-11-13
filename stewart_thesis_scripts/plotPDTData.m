% plotPDTData.m
%  
% From Teo et al (2004): % PAT-depth-temperature (PDT data): 
% During each time interval, the PDT function recorded the min and max 
% temperature experienced at the shallowest adn deepest depths attained 
% by the tagged fish, and at 6 additional depths between those points.  
%  
% Outside Functions Called: 
% plotOxygenProxy.m % J. Stewart
% tagDailyHistos_ShallowPDT % J. Stewart % this is necessary for habitatID_eachTag
%  
% AUTHOR    : J. Stewart jules32@gmail.com 
% DATE      : 21-Mar-2011 11:27:31  
% Revision  : 1.00  
% DEVELOPED : 7.9.0.529 (R2009b) OSX 
% FILENAME  : plotPDTData.m 

Mfilename = mfilename;

%% plot temperature/depth profiles

% temp vs. depth for whole tag
figure;
set(gcf,'Position',[12   213   842   688])
plot(PDTs.TempMean, PDTs.PDTdepth, '.', 'Color', [0 0.5020 0]);
ylabel('Depth [m]', 'FontSize',20, 'FontWeight', 'bold')
xlabel('Temperature [\circC]', 'FontSize',20, 'FontWeight', 'bold')
title({[tagName ' Temperature v. Depth, PDT Data'];[num2str(DeployTot_da,'%2.1f') ' days ']}, 'FontSize',20, 'FontWeight', 'bold')
set(gca, 'YDir', 'reverse');
FigName = [tagNumS 'PDT_tempVdepth.pdf'];
annotate_JS(Mfilename, gcf, FigName)
% orient landscape % save
% print('-dpdf', [dir1 FigName]);

% temp vs. depth daily subplots
Count = 0;
numberofpanels = 6;
numcols = 3;
percx = .6;
percy = .75;
Temptrix = zeros(length(DeployDays), 2);
Depthtrix = zeros(length(DeployDays), 2);

for i = 1:length(DeployDays)
    dlog = DEPLOYINDEX == i;
    datelabel = PDTs.PDTdate(dlog);
    Count = Count + 1;
    fignum = floor((Count - 1)./numberofpanels)+1; % set the number for 6-panel figures
    fignum2 = fignum+100;
    panelnum = rem(Count - 1, numberofpanels)+1; % rem takes the remainder after division
    if panelnum == 1
        figure(fignum2); clf;
    end
    subplot(ceil(numberofpanels/numcols), numcols, panelnum); hold on;
    plot(PDTs.TempMean(dlog), PDTs.PDTdepth(dlog), '.', 'Color', [0 0.5020 0]);
    set(gca,'XLim', [min(PDTs.TempMean) max(PDTs.TempMean)], 'YLim', [min(PDTs.PDTdepth), max(PDTs.PDTdepth)]);
    text(max(PDTs.TempMean)*percx, max(PDTs.PDTdepth)*percy, {['Tmin = ' num2str(min(PDTs.TempMean(dlog))) ' ']; ['Tmax = ' num2str(max(PDTs.TempMean(dlog))) ' ']; ...
        ['Dmin = ' num2str(min(PDTs.PDTdepth(dlog))) ' ']; ['Dmax = ' num2str(max(PDTs.PDTdepth(dlog))) ' ']}, 'fontsize', 12, 'fontweight', 'bold')
    title(datestr(DeployDates(i),6), 'FontSize',12, 'FontWeight', 'bold') 
    set(gca, 'YDir', 'reverse');
    Temptrix(i,:) = [min(PDTs.TempMean(dlog)) max(PDTs.TempMean(dlog))];
    Depthtrix(i,:) = [min(PDTs.PDTdepth(dlog)), max(PDTs.PDTdepth(dlog))];
    if panelnum == 6 
       key = {[tagName ' Daily Temperature/Depth Profiles ' num2str(fignum)]};
       [h3] = suptitle(key);
       [ax,h(1)] = suplabel('Depth [m] ','x'); %('Day of the Year ','x');
       [ax,h(2)] = suplabel('Temperature [\circC] ','y');  
       set(h3,'FontSize',24, 'FontWeight', 'bold')
       set(h,'FontSize',20, 'FontWeight', 'bold')
       FigName = [tagNumS 'PDT_daily_depthVtemp' num2str(fignum) '.pdf'];
      annotate_JS(Mfilename, gcf, FigName)
%        orient landscape % save
%        print('-dpdf', [dir1 FigName]);
   end
end
if panelnum < 6 % catch for when there aren't 6 panels on a figure
    key = {[tagName ' Daily Temperature/Depth Profiles ' num2str(fignum)]};
    [h3] = suptitle(key);
    [ax,h(1)] = suplabel('Temperature [\circC] ','x'); 
    [ax,h(2)] = suplabel('Depth [m] ','y');
    set(h3,'FontSize',24, 'FontWeight', 'bold')
    set(h,'FontSize',20, 'FontWeight', 'bold')
    FigName = [tagNumS 'PDT_daily_depthVtemp' num2str(fignum) '.pdf'];
   annotate_JS(Mfilename, gcf, FigName)
%     orient landscape % save
%     print('-dpdf', [dir1 FigName]);
end

DepthTagMax = Depthtrix(:,2);

%% boxplots as a series

ybit = 800; %dbar(end)+20
    
% depth
figure
subplot(2,1,1)
boxplot(DepthN(~isnan(DepthN)), DEPLOYINDEX(~isnan(DepthN)),'widths', 0.2, 'colors', 'k','symbol', 'ok', 'outliersize', 3);
set(gca,'YLim', [0 ybit], 'YDir', 'reverse')
title('Nighttime','FontSize', 14, 'FontWeight', 'bold')
subplot(2,1,2)
boxplot(DepthD(~isnan(DepthD)), DEPLOYINDEX(~isnan(DepthD)),'widths', 0.2, 'colors', 'k','symbol', 'ok', 'outliersize', 3);
set(gca, 'YLim', [0 ybit], 'YDir', 'reverse') 
title('Daytime','FontSize', 14, 'FontWeight', 'bold')
key = {[tagName ' Daily Depth Boxplots PDT data ' buflab]};
[h3] = suptitle(key);
[ax,h(1)] = suplabel('Day ','x'); %('Day of the Year ','x');
[ax,h(2)] = suplabel('Depth (m)','y');
set(h3,'FontSize',24, 'FontWeight', 'bold')
set(h,'FontSize',20, 'FontWeight', 'bold')
% fixfig
FigName = [tagName '_BoxplotDepthSeriesPDT' buflab '.pdf']; % buflab assigned in processPDTData.m
annotate_JS(Mfilename, gcf, FigName)
% orient landscape % save
% print('-dpdf', [dir1 FigName]);

% temp
figure
subplot(2,1,1)
boxplot(TempN(~isnan(TempN)), DEPLOYINDEX(~isnan(TempN)),'widths', 0.2, 'colors', 'k','symbol', 'ok', 'outliersize', 3);
set(gca,'YLim', [0 tbar(end)+2])
title('Nighttime','FontSize', 14, 'FontWeight', 'bold')
subplot(2,1,2)
boxplot(TempD(~isnan(TempD)), DEPLOYINDEX(~isnan(TempD)),'widths', 0.2, 'colors', 'k','symbol', 'ok', 'outliersize', 3);
set(gca, 'YLim', [0 tbar(end)+2])
title('Daytime','FontSize', 14, 'FontWeight', 'bold')
key = {[tagName ' Daily Temperature Boxplots PDT data ' buflab]};
[h3] = suptitle(key);
[ax,h(1)] = suplabel('Day ','x'); %('Day of the Year ','x');
[ax,h(2)] = suplabel('Temperature (\circC)','y');
set(h3,'FontSize',24, 'FontWeight', 'bold')
set(h,'FontSize',20, 'FontWeight', 'bold')
FigName = [tagName '_BoxplotTempSeriesPDT' buflab '.pdf'];
annotate_JS(Mfilename, gcf, FigName)
% orient landscape % save
% print('-dpdf', [dir1 FigName]);

%% for MBARI comparison

DEPLOYINDEXuni = unique(DEPLOYINDEX);
DEPLOYINDEXday = DEPLOYINDEX(~isnan(DepthD));
DepthDay = DepthD(~isnan(DepthD));
dd = 1:50:2000;

DepthDperc_daily = [];
CCSmaxDtemp = [];
for i = 1:length(DEPLOYINDEXuni)
    ilog = DEPLOYINDEXuni(i) == DEPLOYINDEXday;
    depthtemp = DepthDay(ilog);
    depthhist = hist(depthtemp,dd);
    DepthDperc_daily = [DepthDperc_daily (depthhist./sum(depthhist)*100)']; % DepthDperc_daily is what it is called in processSeriesData.m
    CCSmaxDtemp = [CCSmaxDtemp; max(depthtemp)];
end

% save([CompareTagDir tagName '_DepthHistoDaily.mat'], 'DepthDperc_daily', 'CCSmaxDtemp', 'dd')


%% datetime vs depth with MEAN temp as color

figure(200)
set(gcf,'Position',[12   213   842   688],'Color','w')
%     plot(PDTdate, PDTdepth*-1, ':', 'Color', [.8 .8 .8]); %this is missleading since every hour is split up into 8 bins
%     hold on
scatter(PDTs.PDTdate,PDTs.PDTdepth,20,PDTs.TempMean,'filled');
%     hold off
if exist('DayNight_boxes.m','file') %if the function is present the figure will have night boxes
    [x,y] = DayNight_boxes(PDTs.PDTdate,sr/24,ss/24,get(gca,'Ylim'),0); %will plot boxes
    figure(200)
    set(gcf,'Position',[12   213   842   688],'Color','w')
    fill(x,y,[0.8 0.8 0.8], 'EdgeColor','none')
    hold on
    scatter(PDTs.PDTdate,PDTs.PDTdepth,20,PDTs.TempMean,'filled');
    hold off
end
% set(gca,'XLim',[min(PDT.Datetime) max(PDT.Datetime)])
set(gca, 'YDir', 'reverse')
title([tagName ' PDTdata2'])
ylabel('Depth [m]')
xlabel('Time [GMT]')
Hcbar = colorbar('location','SouthOutside', 'LineWidth', 1);
caxis(cmapscale);%Color axis scaling
xlabel(Hcbar,'Temperature [\circC]');
datetick('x',6,'keeplimits')
FigName = [tagNumS 'PDT_timeVdepth_tempMean.pdf'];
annotate_JS(Mfilename, gcf, FigName);
% orient landscape % save
% print('-dpdf', [dir1 FigName]);

%% datetime vs depth with MIN,MAX temp as color

figure(200)
set(gcf,'Position',[12   213   842   688],'Color','w')
%     plot(PDTdate, PDTdepth*-1, ':', 'Color', [.8 .8 .8]); %this is missleading since every hour is split up into 8 bins
%     hold on
scatter(PDTs.PDTdate,PDTs.PDTdepth,20,PDTs.TempMean,'filled');
%     hold off
if exist('DayNight_boxes.m','file') %if the function is present the figure will have night boxes
    [x,y] = DayNight_boxes(PDTs.PDTdate,sr/24,ss/24,get(gca,'Ylim'),0); %will plot boxes
    figure(200)
    set(gcf,'Position',[12   213   842   688],'Color','w')
    fill(x,y,[0.8 0.8 0.8], 'EdgeColor','none')
    hold on
    scatter(PDTs.PDTdate,PDTs.PDTdepth,20,PDTs.TempMin,'filled');
    hold off
end
% set(gca,'XLim',[min(PDT.Datetime) max(PDT.Datetime)])
set(gca, 'YDir', 'reverse')
title([tagName ' PDTdata2'])
ylabel('Depth [m]')
xlabel('Time [GMT]')
Hcbar = colorbar('location','SouthOutside', 'LineWidth', 1);
caxis(cmapscale);%Color axis scaling
xlabel(Hcbar,'Temperature [\circC]');
datetick('x',6,'keeplimits')
FigName = [tagNumS '_TimeSeries_TimeDepth_TempMean.pdf'];
annotate_JS(Mfilename, gcf, FigName);
% orient landscape % save
% print('-dpdf', [dir1 FigName]);


%%

tagDailyHistos_ShallowPDT % J. Stewart % this is necessary for habitatID_eachTag
% plotOxygenProxy % J. Stewart 

%% 
  
disp('Completed plotPDTData.m') 
% ===== EOF [plotPDTData.m] ======  
