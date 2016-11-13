% compareTagsFigs.m
%
% Called From: compareTags.m

% Description: makes all figures associated with compareTags
%
% Outside Functions Called:
% sunrise_.m M. Broenkow
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 07-Oct-2011 15:40:54
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : compareTagsFigs.m

Mfilename = mfilename;
clear H % otherwise legend doesn't work
cd(dirCompare);

%% Make min/max plots for depth (blue) and temp (green)

Mfilename = mfilename;
DepthMaxes = zeros(length(TagIDs),1);

% temperature figure
figure
for i = 1:length(TagIDs)
    tlog = CompareDT(:,7) == TagIDs(i);
    TagIDtempA = CompareDT(tlog,2);
    TagIDtemp = TagIDtempA(1);
    
    DepthMaxes(i) = max(CompareDT(tlog,6));
    
    subplot(rows,cols,i)
    h(2) = area(CompareDT(tlog,1), CompareDT(tlog,4)); % max temp
    hold on
    h(1) = area(CompareDT(tlog,1), CompareDT(tlog,3)); % min temp
    set(gca, 'YLim', [0 max(CompareDT(:,4))]);
    set(h(2),'FaceColor',[0.76 0.87 0.78])
    set(h(1),'FaceColor',[1 1 1])
    set(gca,'Layer','top')
    plot(CompareDT(tlog,1), CompareDT(tlog,3:4), 'k','LineWidth',2)
    set(gca, 'fontsize', 12, 'fontweight', 'bold')
    title({[num2str(floor(TagIDtemp)) '-' datestr(min(CompareDT(tlog,1)),11) ': ' datestr(min(CompareDT(tlog,1)),6) '-' datestr(max(CompareDT(tlog,1)),6)];...
        [num2str(sum(tlog)) ' days ']}, 'FontSize',12, 'FontWeight', 'bold'); %[num2str(sum(tlog)) ' days, ' num2str(Rates_kmday(i)) ' km/day']}, 'FontSize',12, 'FontWeight', 'bold')
    set(gca,'XTick', CompareDT(tlog,1))
    datetick('x','dd','keepticks')
end
key = {'Daily Min/Max Temperature Tag Comparisons (Time Series) '};
% [h3] = suptitle(key);
% [ax,h(1)] = suplabel('Date ','x');
% [ax,h(2)] = suplabel('Temperature (\circC) ','y');
% set(h3,'FontSize',24, 'FontWeight', 'bold')
% set(h,'FontSize',20, 'FontWeight', 'bold')
% FigName = ['SeriesCompareMinMaxT_' locID '.pdf'];
% annotate_JS(Mfilename, gcf, FigName)
% orient landscape % save
% print('-dpdf', FigName);

% depth figure
figure
for i = 1:length(TagIDs)
    tlog = CompareDT(:,7) == TagIDs(i);
    TagIDtempA = CompareDT(tlog,2);
    TagIDtemp = TagIDtempA(1);
    subplot(rows,cols,i)
    h(2) = area(CompareDT(tlog,1), CompareDT(tlog,6)); % max depth
    hold on
    h(1) = area(CompareDT(tlog,1), CompareDT(tlog,5)); % min depth
    set(gca, 'YLim', [0 max(CompareDT(:,6))]);
    set(gca, 'YDir', 'reverse')
    set(h(2),'FaceColor',[0.7 0.78 1])
    set(h(1),'FaceColor',[1 1 1])
    set(gca,'Layer','top')
    plot(CompareDT(tlog,1), CompareDT(tlog,5:6), 'k','LineWidth',2)
    set(gca, 'fontsize', 12, 'fontweight', 'bold')
    title({[num2str(TagIDtemp) '-' datestr(min(CompareDT(tlog,1)),11) ': ' datestr(min(CompareDT(tlog,1)),6) '-' datestr(max(CompareDT(tlog,1)),6)],...
        [num2str(sum(tlog)) ' days ']}, 'FontSize',12, 'FontWeight', 'bold') %[num2str(sum(tlog)) ' days, ' num2str(Rates_kmday(i)) ' km/day']}, 'FontSize',12, 'FontWeight', 'bold')
    set(gca,'XTick', CompareDT(tlog,1))
    datetick('x','dd','keepticks')
end
key = {'Daily Min/Max Depth Tag Comparisons (Time Series) '};
% [h3] = suptitle(key);
% [ax,h(1)] = suplabel('Date ','x');
% [ax,h(2)] = suplabel('Depth (m) ','y');
% set(h3,'FontSize',24, 'FontWeight', 'bold')
% set(h,'FontSize',20, 'FontWeight', 'bold')
% FigName = ['SeriesCompareMinMaxD_' locID '.pdf'];
% annotate_JS(Mfilename, gcf, FigName)
% orient landscape % save
% print('-dpdf', FigName);

clear ax

%% make Depth and Temp Boxplots and Histograms and Tables!

for p = 1%:2
    if p == 1 % Depth
        Data = BoxPlotDepthtrix;
        DataLabel = 'Depth-m';
        ydir = 'reverse';
        Bar= dbar;
        DataXLim = 1500;
        DataYLim = 40;
    elseif p == 2 % Temp
        Data = BoxPlotTemptrix;
        DataLabel = 'Temperature-C';
        ydir = 'normal';
        Bar=tbar;
        DataXLim = 30;
        DataYLim = 30;
    end
    
    Bin = Bar(2) - Bar(1);
    barend = max(Bar); 
    dbarAll = 0:Bin:max(Bar);
    
    TagIDstats = []; % for the table below
    TagAllstats = [];
    
    % setup for boxplots
    DataDay = Data(:,4);
    DataNight = Data(:,3);
    DatamaxD = max(DataDay);
    DatamaxN = max(DataNight);
    
    figure
    subplot(1,2,1)
    boxplot(DataNight, Data(:,5), 'colors', 'k','symbol', '.k', 'outliersize', 5); % plot 5 tags
    %     boxplot([DataNight; DataNight], [Data(:,5); repmat(7,length(DataNight),1)], 'colors', 'k','symbol', '.k', 'outliersize', 3); % plot individs and all tags
    set(gca, 'YLim', [0-Bin DataXLim], 'YDir', ydir) % since boxplots are still on y axis
    title('Nighttime','FontSize', 14, 'FontWeight', 'bold')
    subplot(1,2,2)
    boxplot(DataDay, Data(:,5), 'colors', 'k', 'symbol', '.k', 'outliersize', 5); % plot 5 tags
    %     boxplot([DataDay; DataDay], [Data(:,5); repmat(7,length(DataDay),1)], 'colors', 'k','symbol', '.k', 'outliersize', 3); % plot individs and all tags
    set(gca, 'YLim', [0-Bin DataXLim],'YDir', ydir)
    title('Daytime','FontSize', 14, 'FontWeight', 'bold')
    key = [DataLabel ' Boxplots ' locID vertID];
%     [h3] = suptitle(key);
%     [ax,h(1)] = suplabel('Tag Number ','x'); %('Day of the Year ','x');
%     [ax,h(2)] = suplabel(DataLabel, 'y');
%     set(h3,'FontSize',24, 'FontWeight', 'bold')
%     set(h,'FontSize',20, 'FontWeight', 'bold')
%     FigName = ['BoxplotSubplot_' DataLabel '_' locID vertID '.pdf'];
%     annotate_JS(Mfilename, gcf, FigName)
%     fixfig
%     orient landscape % save
%     print('-dpdf', [FigName]);
    
    %make summary information matrices for day and night information.
    SummtrixDay = [];
    SummtrixNight = [];
    nanMedians = [];
    nanStds = [];
    
    for i = 1:length(TagIDs)
        tlog = Data(:,5) == TagIDs(i); % for some reason TagIDs has extra decimal points so it won't be recognized..
        ddtemp = DataDay(tlog);
        summday = summary(ddtemp);
        medtemp = nanmedian(DataDay(tlog));
        summnight = summary(DataNight(tlog));
        nMtemp = nanmedian(DataDay(tlog));
        nanstdday = nanstd(DataDay(tlog));
        nanstdnight = nanstd(DataNight(tlog));
        SummtrixDay = [SummtrixDay summday];
        SummtrixNight = [SummtrixNight summnight];
        nanMedians = [nanMedians nMtemp];
        nanStds = [nanStds; nanstdday nanstdnight];
    end
    nanSDs = [nanStds(:,2)' nanStds(:,1)'];
    
%     %Stat output % commented out 11/7/13
%     %Daytime
%     SummtrixDay(4,:)' % median range
%     nanmedian(DataDay) %median of all tags
%     nanmean(DataDay)
%     nanstd(DataDay)
%     iqr(DataDay(~isnan(DataDay)))
%     
%     %Nighttime
%     SummtrixNight(4,:)' % median range
%     nanmedian(DataNight) %median of all tags
%     nanmean(DataNight)
%     nanstd(DataNight)
%     iqr(DataNight(~isnan(DataNight)))
    
    % combined histograms for depth and temperature
    % setup
    dmindex = 1;
    dmaxdex = DataYLim/Bin + 1;
    TimeinSwath = [];
    Data_binnedN = hist(DataDay,dbarAll);%count within each bin
    percDNalltags = ((Data_binnedN/sum(Data_binnedN))*100)';
    Data_binnedD = hist(DataNight,dbarAll); %count within each bin
    percDDalltags = ((Data_binnedD/sum(Data_binnedD))*100)';
    
%     summdepthday = summary(DataDay);
%     summdepthnight = summary(DataNight);
%     summdepthday(5)-summdepthday(3)
%     
    %%%% histograms%%%%%
    figure;
    H(2) = bar(dbarAll,percDDalltags,'BarWidth',1);
    hold on
    H(1) = bar(dbarAll,percDNalltags,'BarWidth',1);
    hold off
    set(H(1),'facecolor','none','EdgeColor','k','LineWidth',1.5) %night
    set(H(2),'facecolor',[0.5 0.5 0.5],'EdgeColor',[0.2 0.2 0.2],'LineWidth',1.5) %day
    set(gca, 'XLim', [0-Bin DataXLim], 'YLim', [0 DataYLim])
    ylabel('Percent of Time','FontSize', 14, 'FontWeight', 'bold')
    xlabel(DataLabel,'FontSize', 14, 'FontWeight', 'bold')
    legend(H,'Day','Night','Location','NorthEast')
    title({[locID vertID ' Combined Tag Histograms: ' DataLabel]; ['Daytime mean = ' num2str(summdepthday(1))...
        ', sd = ' num2str(nanstd(DataDay)) ', median = ' num2str(summdepthday(4))];...
        ['Nighttime mean = ' num2str(summdepthnight(1)) ', sd = ' num2str(nanstd(DataNight)) ...
        ', median = ' num2str(summdepthnight(4))]});
%     FigName = 'TagsDaytimeDepth_CA_4tags2009.pdf';%['HistosAllTags_' DataLabel '_' locID vertID '.pdf']; % when run with just forage tags, label it.
%     annotate_JS(Mfilename, gcf, FigName)
%     %     orient landscape % save
    %     print('-dpdf', [FigName]);
    
    % make for fig in c_ThesisFigs for ch3. Run compareTags.m with only 4 2009
    % archival tags (#3-6), and above in this section set p=1. and CHANGE dbar in compareTags to :50:
%     CA_4tagDayBins = dbarAll;
%     CA_4tagDayPerc = percDNalltags;
%     CA_4tagMaxTempDepth = CompareDT(:,1:7);
%     compareTagsMBARI

% make for fig in c_ThesisFigs for ch3. Run compareTags.m with only 4 2009+1 2008
%     CA_5tagDayBins = dbarAll;
%     CA_5tagDayPerc = percDNalltags;
%     CA_5tagMaxTempDepth = CompareDT(:,1:7);
%     save -mat '/Users/juliastewart/Dropbox/Thesis/Dg_Abundance/MBARI_Dg/CA_5tagsDayHist.mat'...
%     CA_5tagDayBins CA_5tagDayPerc CA_5tagMaxTempDepth...
    
    %     figure;
    %     subplot(1,2,1)
    %     H(2) = bar(dbarAll,percDDalltags,'BarWidth',1);
    %     hold on
    %     H(1) = bar(dbarAll,percDNalltags,'BarWidth',1);
    %     hold off
    %     set(H(1),'facecolor','none','EdgeColor','k','LineWidth',1.5) %night
    %     set(H(2),'facecolor',[0.5020    0.5020    0.5020],'EdgeColor',[0.2 0.2 0.2],'LineWidth',1.5) %day
    %     FigName = ['HistosAllTags_' DataLabel '_' locID vertID '.pdf']; % when run with just forage tags, label it.
    %     annotate_JS(Mfilename, gcf, FigName)
    %     set(gca, 'XLim', [0 2000], 'YLim', [0 25]);
    %     ylabel('CA Tag % at Depth')
    %     FigName = ['CAtags_toCompareVars.pdf']; % when run with just forage tags, label it.
    %     annotate_JS(Mfilename, gcf, FigName)
    %     orient landscape % save
    %     print('-dpdf', [FigName]);
    
    
    
    
    %%%%%TABLES%%%%%%
    [ColNamesDay] = makeLabels(TagIDs, ' day'); % A. Booth
    [ColNamesNight] = makeLabels(TagIDs, ' night'); % A. Booth
    ColNames = ['stats'; ColNamesNight; 'all tags'; ColNamesDay; 'all tags']';
    RowNames = {'mean', 'min', '1st quarter', 'median', '3rd quarter', 'max', 'sd'}';
    
    TagIDstats = [TagIDstats SummtrixNight SummtrixDay];
    TagAllstats = [TagAllstats summdepthnight summdepthday];
    TagIDstatsM = [TagIDstats(1:end-1,:); nanSDs];
    TagAllstatsM = [TagAllstats(1:end-1,:); nanstd(DataNight) nanstd(DataDay)];
    
    figure
    tbl = table([ColNames], RowNames, ...
        ceil([TagIDstatsM(:,1:length(TagIDs)) TagAllstatsM(:,1)...
        TagIDstatsM(:,1+length(TagIDs):length(TagIDs)+length(TagIDs)) TagAllstatsM(:,2)])); % fix for now to GOC data can fit
    plot_table(tbl);
    [x,y] = getAxisInset( .1, .95 );
    pos = [360   502   670   420];
    set( gcf,'pos', pos );
    text(x,y, [locID vertID ' Tags ' DataLabel ' Summaries '], 'fontsize', 18,'fontweight', 'bold');
%     FigName = ['Table_' DataLabel '_' locID vertID '.pdf']; % when run with just forage tags, label it.
%     %     fixfig
%     annotate_JS(Mfilename, gcf, FigName)
%     orient landscape % save
%     print('-dpdf', [FigName]);
end

%%

disp('Completed compareTagsFigs.m')
% ===== EOF [compareTagsFigs.m] ======
