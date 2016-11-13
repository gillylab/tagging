% plotSeriesData.m
% plots a lot, and also makes a structure of Profiles for later.
%
% called from processSeriesData.m
%
% Outside Functions Called:
% ps2pdf.m % Matlab Central
% plotSeriesDataPublication.m % J. Stewart
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 24-Mar-2011 16:29:28
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : plotSeriesData.m

Mfilename = mfilename;

%% Plot Series Histograms

summdepthday = summary(DepthTag(isfinite(timeD)));
summdepthnight = summary(DepthTag(isfinite(timeN)));
summtempday = summary(TempTag(isfinite(timeD)));
summtempnight = summary(TempTag(isfinite(timeN)));

DepthHistotrix = [repmat(tagNumD,length(percDN),1) repmat(yr,length(percDN),1) dbar' percDN percDD];
TempHistotrix = [repmat(tagNumD,length(percTN),1) repmat(yr,length(percTN),1) tbar' percTN percTD];

if pf
    % depth
    figure;
    H(2) = barh(dbar,percDD,'BarWidth',1);
    hold on
    H(1) = barh(dbar,percDN,'BarWidth',1);
    hold off
    set(H(1),'facecolor','none','EdgeColor','k','LineWidth',1.5) %night
    set(H(2),'facecolor',[0.5020    0.5020    0.5020],'EdgeColor',[0.2 0.2 0.2],'LineWidth',1.5) %dayxlabel('Percent of Time','FontSize', 14, 'FontWeight', 'bold')
    ylabel(['Depth Bin (' num2str(binsizeD) 'm)'],'FontSize', 14, 'FontWeight', 'bold')
    legend(H,'Night','Day','Location','SouthEast')
    set(gca,'YLim',[0,dbar(end)+20], 'YDir', 'reverse')%[min(dbar), max(dbar)])
    set(gca, 'fontsize', 12, 'fontweight', 'bold')
    box off
    h1 = gca;
    title({[tagName ' Depths, ' num2str(DeployTot_da,'%2.1f') ' days']; ['Day median = ' num2str(summdepthday(4)) ', Night median = ' num2str(summdepthnight(4))];...
        [datestr(min(DeployDates),22) ' - ' datestr(max(DeployDates),22) ' (GMT)']}, 'fontsize', 16, 'fontweight', 'bold');
    FigName = [tagName '_HistoDepth_bins' num2str(binsizeD) 'm.pdf'];
    annotate_JS(Mfilename, gcf, FigName);
%     orient landscape % save
%     print('-dpdf', [dir1 FigName]);
    
    % temperature
    figure;
    H(2) = barh(tbar,percTD,'BarWidth',1);
    hold on
    H(1) = barh(tbar,percTN,'BarWidth',1);
    hold off
    set(H(1),'facecolor','none','EdgeColor','k','LineWidth',1.5) %night
    set(H(2),'facecolor',[0.5020    0.5020    0.5020],'EdgeColor',[0.2 0.2 0.2],'LineWidth',1.5) %day
    xlabel('Percent of Time','FontSize', 14, 'FontWeight', 'bold')
    ylabel(['Temperature Bin (' num2str(binsizeT) '\circC)'],'FontSize', 14, 'FontWeight', 'bold')
    legend(H,'Night','Day','Location','SouthEast')
    set(gca,'YLim',[0,tbar(end)+1])%[min(dbar), max(dbar)])
    set(gca, 'fontsize', 12, 'fontweight', 'bold')
    box off
    h1 = gca;
    title({[tagName ' Temperatures, ' num2str(DeployTot_da,'%2.1f') ' days']; ['Day median = ' num2str(summtempday(4)) ', Night median = ' num2str(summtempnight(4))];...
        ['2.5 percentile = ' num2str(prctile(TempTag(isfinite(timeD)), 2.5)) ', 97.5 percentile = ' num2str(prctile(TempTag(isfinite(timeD)), 97.5))];...
        [datestr(min(DeployDates),22) ' - ' datestr(max(DeployDates),22) ' (GMT)']}, 'fontsize', 16, 'fontweight', 'bold');
    FigName = [tagName '_HistoTemp_bins' num2str(binsizeT) 'C.pdf'];
    annotate_JS(Mfilename, gcf, FigName);
%     orient landscape % save
%     print('-dpdf', [dir1 FigName]);
    
    %% Depth histos with OML--should this just go to Tags_CTD?
    
    if OMLhistfig %defined in importTagData_csv
        figure
        H(2) = barh(dbar,percDD,'BarWidth',1);
        hold on
        H(1) = barh(dbar,percDN,'BarWidth',1);
        set(H(1),'facecolor','none','EdgeColor','k','LineWidth',1.5) %night
        set(H(2),'facecolor',[0.5020    0.5020    0.5020],'EdgeColor',[0.2 0.2 0.2],'LineWidth',1.5) %day
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
    for i = 1:length(DeployDays)
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
        set(H(1),'facecolor','none','EdgeColor','k','LineWidth',1.5) %night
        set(H(2),'facecolor',[0.5020    0.5020    0.5020],'EdgeColor',[0.2 0.2 0.2],'LineWidth',1.5) %daytitle(['Day ' num2str(i) ', ' datestr(DeployDates(i),6) ', ' num2str(ceil(DeployHours(i))) ' hours'], 'FontSize',12, 'FontWeight', 'bold')
        legend(H,'Day','Night','Location','SouthEast')
        set(gca,'YLim',[0, MaxDepthLim], 'YDir', 'reverse', 'XLim', [0 100]);
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
            %         orient landscape % save
            %         print('-dpdf', [dir1 FigName]);
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
        %     orient landscape % save
        %     print('-dpdf', [dir1 FigName]);
    end
    
    %% more histos
    
    % temp
    Count = 0;
    for i = 1:length(DeployDays)
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
        set(H(1),'facecolor','none','EdgeColor','k','LineWidth',1.5) %night
        set(H(2),'facecolor',[0.5020    0.5020    0.5020],'EdgeColor',[0.2 0.2 0.2],'LineWidth',1.5) %daytitle(['Day ' num2str(i) ', ' datestr(DeployDates(i),6) ', ' num2str(ceil(DeployHours(i))) ' hours'], 'FontSize',12, 'FontWeight', 'bold')
        legend(H,'Day','Night','Location','SouthEast')
        set(gca,'YLim',[0,tbar(end)+1], 'XLim', [0 100]);
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
            %         orient landscape % save
            %         print('-dpdf', [dir1 FigName]);
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
        %     orient landscape % save
        %     print('-dpdf', [dir1 FigName]);
    end
    
    %% daily boxplots
    
    % depth
    Count = 0;
    for i = 1:length(DeployDays)
        dlog = DEPLOYINDEX == DeployDays(i); %i seems to work for GOC, the other way was fine with CA. It's when it spans a month...DeployDays(i);DeployDays(i);
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
            key = {[tagName ' Daily Depth Boxplots ' buflab '_' num2str(fignum)]};
            [h3] = suptitle(key);
            [ax,h] = suplabel('Depth (m)','y');
            set(h3,'FontSize',24, 'FontWeight', 'bold')
            set(h,'FontSize',20, 'FontWeight', 'bold')
            FigName = [tagName '_BoxplotDepth_' buflab '_Daily_' num2str(fignum) '.pdf'];
            annotate_JS(Mfilename, gcf, FigName)
%             orient landscape % save
%             print('-dpdf', [dir1 FigName]);
        end
    end
    if panelnum < 6 % catch for when there aren't 6 panels on a figure
        %     fixfig
        key = {[tagName ' Daily Depth Boxplots ' buflab '_' num2str(fignum)]};
        [h3] = suptitle(key);
        [ax,h] = suplabel('Depth (m)','y');
        set(h3,'FontSize',24, 'FontWeight', 'bold')
        set(h,'FontSize',20, 'FontWeight', 'bold')
        FigName = [tagName '_BoxplotDepth_' buflab '_Daily_' num2str(fignum) '.pdf'];
        annotate_JS(Mfilename, gcf, FigName)
%         orient landscape % save
%         print('-dpdf', [dir1 FigName]);
    end
    
    % temp
    Count = 0;
    for i = 1:length(DeployDays)
        dlog = DEPLOYINDEX == DeployDays(i); % i
        TempNn = TempN(dlog);
        TempDd = TempD(dlog);
        if length(TempNn) > length(TempDd) % vectors must be the same length
            TempDD = [TempDd; nan((length(TempNn)-length(TempDd)),1)];
            TempNN = TempNn;
        else
            TempDD = TempDd;
            TempNN = [TempNn; nan((length(TempDd)-length(TempNn)),1)];
        end
        boxplotND = [TempNN TempDD];
        Count = Count + 1;
        fignum = floor((Count - 1)./numberofpanels)+1; % set the number for 6-panel figures
        fignum2 = fignum+260;
        panelnum = rem(Count - 1, numberofpanels)+1; % rem takes the remainder after division
        if panelnum == 1
            figure(fignum2); clf;
        end
        subplot(ceil(numberofpanels/numcols), numcols, panelnum); hold on;
        boxplot(boxplotND, 'colors', 'k', 'symbol', '.k', 'outliersize', 5, 'labels', {'Night', 'Day'});
        set(gca, 'YLim', [0 tbar(end)+2])
        title([datestr(DeployDates(i),6) ', ' num2str(ceil(DeployHours(i))) ' hours'], 'FontSize',12, 'FontWeight', 'bold')
        if panelnum == 6
            %         fixfig
            key = {[tagName ' Daily Temperature Boxplots ' buflab '_' num2str(fignum)]};
            [h3] = suptitle(key);
            [ax,h] = suplabel('Temperature (\circC)','y');
            set(h3,'FontSize',24, 'FontWeight', 'bold')
            set(h,'FontSize',20, 'FontWeight', 'bold')
            FigName = [tagName '_BoxplotTemp_' buflab '_Daily' num2str(fignum) '.pdf'];
            annotate_JS(Mfilename, gcf, FigName)
%             orient landscape % save
%             print('-dpdf', [dir1 FigName]);
        end
    end
    if panelnum < 6 % catch for when there aren't 6 panels on a figure
        %     fixfig
        key = {[tagName ' Daily Temperature Boxplots ' buflab '_' num2str(fignum)]};
        [h3] = suptitle(key);
        [ax,h] = suplabel('Temperature (\circC)','y');
        set(h3,'FontSize',24, 'FontWeight', 'bold')
        set(h,'FontSize',20, 'FontWeight', 'bold')
        FigName = [tagName '_BoxplotTemp_' buflab '_Daily' num2str(fignum) '.pdf'];
        annotate_JS(Mfilename, gcf, FigName)
%         orient landscape % save
%         print('-dpdf', [dir1 FigName]);
    end
    
    %% boxplots as a series
    
    % depth
    figure
    subplot(2,1,1)
    boxplot(DepthN(~isnan(DepthN)), DEPLOYINDEX(~isnan(DepthN))-(DeployDays(1)-1),'widths', 0.3, 'colors', 'k','symbol', 'ok', 'outliersize', 3);
    set(gca,'YLim', [0 dbar(end)+20], 'YDir', 'reverse')
    title('Nighttime','FontSize', 14, 'FontWeight', 'bold')
    subplot(2,1,2)
    boxplot(DepthD(~isnan(DepthD)), DEPLOYINDEX(~isnan(DepthD))-(DeployDays(1)-1),'widths', 0.3, 'colors', 'k','symbol', 'ok', 'outliersize', 3);
    set(gca, 'YLim', [0 dbar(end)+20], 'YDir', 'reverse')
    title('Daytime','FontSize', 14, 'FontWeight', 'bold')
    key = {[tagName ' Daily Depth Boxplots ' buflab]};
    [h3] = suptitle(key);
    [ax,h(1)] = suplabel('Day ','x'); %('Day of the Year ','x');
    [ax,h(2)] = suplabel('Depth (m)','y');
    set(h3,'FontSize',24, 'FontWeight', 'bold')
    set(h,'FontSize',20, 'FontWeight', 'bold')
    % fixfig
    FigName = [tagName '_BoxplotDepthSeries' buflab '.pdf']; % buflab assigned in processSeriesData.m
    annotate_JS(Mfilename, gcf, FigName)
%     orient landscape % save
%     print('-dpdf', [dir1 FigName]);
    
    % temp
    figure
    subplot(2,1,1)
    boxplot(TempN(~isnan(TempN)), DEPLOYINDEX(~isnan(TempN))-(DeployDays(1)-1),'widths', 0.3, 'colors', 'k','symbol', 'ok', 'outliersize', 3);
    set(gca,'YLim', [0 tbar(end)+2])
    title('Nighttime','FontSize', 14, 'FontWeight', 'bold')
    subplot(2,1,2)
    boxplot(TempD(~isnan(TempD)), DEPLOYINDEX(~isnan(TempD))-(DeployDays(1)-1),'widths', 0.3, 'colors', 'k','symbol', 'ok', 'outliersize', 3);
    set(gca, 'YLim', [0 tbar(end)+2])
    title('Daytime','FontSize', 14, 'FontWeight', 'bold')
    key = {[tagName ' Daily Temperature Boxplots ' buflab]};
    [h3] = suptitle(key);
    [ax,h(1)] = suplabel('Day ','x'); %('Day of the Year ','x');
    [ax,h(2)] = suplabel('Temperature (\circC)','y');
    set(h3,'FontSize',24, 'FontWeight', 'bold')
    set(h,'FontSize',20, 'FontWeight', 'bold')
    FigName = [tagName '_BoxplotTempSeries' buflab '.pdf'];
    annotate_JS(Mfilename, gcf, FigName)
%     orient landscape % save
%     print('-dpdf', [dir1 FigName]);
    
end

%% For NSF 2013 report (just above where this crashes because of JS_annotate)

DN_melt = reshape(DepthNbinned_daily, [],1);
DD_melt = reshape(DepthDbinned_daily, [],1);
TN_melt = reshape(TempNbinned_daily, [],1);
TD_melt = reshape(TempDbinned_daily, [],1);

dbar_melt = repmat(dbar', size(DepthNbinned_daily, 2),1);
tbar_melt = repmat(tbar', size(TempNbinned_daily, 2),1);

% put deployment day indicator
D2trx = [];
for p = 1:size(DepthNbinned_daily, 2)
    p2 = repmat(p, size(DepthNbinned_daily, 1),1);
    D2trx = [D2trx; p2];
end

T2trx = [];
for p = 1:size(TempNbinned_daily, 2)
    p2 = repmat(p, size(TempNbinned_daily, 1),1);
    T2trx = [T2trx; p2];
end

Dtocsv = [dbar_melt, D2trx, DN_melt, DD_melt];
Ttocsv = [tbar_melt, T2trx, TN_melt, TD_melt];

Hs = cellstr(['depth_m  '; 'deployday'; 'night    '; 'day      ']);

if TagLocation == 0 
    id = 'GOC-'; 
elseif TagLocation == 1
    id = 'CCS-';
end

csvwrite_with_headers(['/Users/jstewart/Dropbox/Thesis/Dg_Tagging/_NSF_report2013/' id tagName '_counts_Depth_50m_bins.csv'], Dtocsv, Hs, 0,0)
csvwrite_with_headers(['/Users/jstewart/Dropbox/Thesis/Dg_Tagging/_NSF_report2013/' id tagName '_counts_Temp_point5C_bins.csv'], Ttocsv, Hs, 0,0)
disp('')
disp('¡¡¡ GOC SquidNumbers are different in DSR paper!!!!!!!!!! Will need to redo by hand')


%% temp vs. depth. also save as structure.

% temp vs. depth for whole tag
figure;
set(gcf,'Position',[12   213   842   688])
plot(TempTag, DepthTag, '.', 'Color', [0 0.5020 0]);
ylabel('Depth [m]', 'FontSize',20, 'FontWeight', 'bold')
xlabel('Temperature [\circC]', 'FontSize',20, 'FontWeight', 'bold')
title([tagName ' Temperature v. Depth, SeriesData'], 'FontSize',20, 'FontWeight', 'bold')
set(gca, 'YDir', 'reverse');
FigName = [tagName '_Profile_TempDepth.pdf'];
annotate_JS(Mfilename, gcf, FigName)
% orient landscape % save
% print('-dpdf', FigName);

% temp vs. depth daily subplots
Count = 0;
percx = .6;
percy = .75;

Temptrix = zeros(length(DeployDays), 2); % just the min and max temp and depths
Depthtrix = zeros(length(DeployDays), 2);

for i = 1:length(DeployDays)
    dlog = DEPLOYINDEX == DeployDays(i); %==i seems to work for BLUE GOC, ==DeployDays(i) works with CA. It's when it spans a month...
    Count = Count + 1;
    fignum = floor((Count - 1)./numberofpanels)+1; % set the number for 6-panel figures
    fignum2 = fignum+100;
    panelnum = rem(Count - 1, numberofpanels)+1; % rem takes the remainder after division
    if panelnum == 1
        figure(fignum2); clf;
    end
    subplot(ceil(numberofpanels/numcols), numcols, panelnum); hold on;
    plot(TempD(dlog), DepthD(dlog), '.', 'Color', [0 0.5020 0]);
    plot(TempN(dlog), DepthN(dlog), '.', 'Color', [0 0.5020 0]);
    set(gca,'XLim', [min(TempTag) max(TempTag)], 'YLim', [min(DepthTag), max(DepthTag)]);
    %     text(max(TempTag)*percx, max(DepthTag)*percy, {['Tmin = ' num2str(min(TempTag(dlog))) ' ']; ['Tmax = ' num2str(max(TempTag(dlog))) ' ']; ...
    %         ['Dmin = ' num2str(min(DepthTag(dlog))) ' ']; ['Dmax = ' num2str(max(DepthTag(dlog))) ' ']}, 'fontsize', 12, 'fontweight', 'bold')
    title([datestr(DeployDates(i),6) ', ' num2str(ceil(DeployHours(i))) ' hours'], 'FontSize',12, 'FontWeight', 'bold')
    set(gca, 'YDir', 'reverse');
    if panelnum == 6
        key = {[tagName ' Daily Temperature/Depth Profiles ' num2str(fignum)]};
        [h3] = suptitle(key);
        [ax,h(1)] = suplabel('Depth [m] ','x'); %('Day of the Year ','x');
        [ax,h(2)] = suplabel('Temperature [\circC] ','y');
        set(h3,'FontSize',24, 'FontWeight', 'bold')
        set(h,'FontSize',20, 'FontWeight', 'bold')
        FigName = [tagName '_Profile_TempDepth_Daily' num2str(fignum) '.pdf'];
        annotate_JS(Mfilename, gcf, FigName)
        %         orient landscape % save
        %         print('-dpdf', FigName);
    end
    
    % save things for later
    Temptrix(i,:) = [min(TempTag(dlog)) max(TempTag(dlog))];
    Depthtrix(i,:) = [min(DepthTag(dlog)), max(DepthTag(dlog))];
    SquidProfile(i).tempAscent = TempAscent(dlog);
    SquidProfile(i).tempDescent = TempDescent(dlog);
    SquidProfile(i).depthAscent = DepthAscent(dlog);
    SquidProfile(i).depthDescent = DepthDescent(dlog);
    %     SquidProfile(i).depthA_Thermocline10 = DepthAscent(TempAscent(dlog)>=9.5 & TempAscent(dlog)<10.5); % for thermocline
    %     SquidProfile(i).depthD_Thermocline10 = DepthDescent(TempDescent(dlog)>=9.5 & TempDescent(dlog)<10.5);
    SquidProfile(i).depthA_Thermocline10 = DepthD(TempD(dlog)==10);%>=9.5 & TempD(dlog)<10.5); % for thermocline
    SquidProfile(i).depthD_Thermocline10 = DepthN(TempN(dlog)==10);%>=9.5 & TempN(dlog)<10.5);
    SquidProfile(i).tempDay = TempD(dlog);
    SquidProfile(i).tempNight = TempN(dlog);
    SquidProfile(i).depthDay = DepthD(dlog);
    SquidProfile(i).depthNight = DepthN(dlog);
    SquidProfile(i).deploydate = DeployDates(i);
    SquidProfile(i).deployhours = DeployHours(i);
end
if panelnum < 6 % catch for when there aren't 6 panels on a figure
    key = {[tagName ' Daily Temperature/Depth Profiles ' num2str(fignum)]};
    [h3] = suptitle(key);
    [ax,h(1)] = suplabel('Temperature [\circC] ','x');
    [ax,h(2)] = suplabel('Depth [m] ','y');
    set(h3,'FontSize',24, 'FontWeight', 'bold')
    set(h,'FontSize',20, 'FontWeight', 'bold')
    FigName = [tagName '_Profile_TempDepth_Daily' num2str(fignum) '.pdf'];
    annotate_JS(Mfilename, gcf, FigName)
    %     orient landscape % save
    %     print('-dpdf', FigName);
end

% filenameSave = ['SquidProfileStructure_' tagNumS '.mat'];
% eval (['save -mat ' filenameSave ' SquidProfile'])

%% test plot thermocline

% figure; hold on
% for i = 1:length(DeployDays)
%     plotmeanA = mean(SquidProfile(i).depthA_Thermocline10);
%     plotmeanD = mean(SquidProfile(i).depthD_Thermocline10);
%     plot(i,plotmeanA, '.k', i, plotmeanD, '.m')
%     set(gca, 'YDir', 'reverse');
% end
%
% for i = 1:length(DeployDays)
%     plot(repmat(i, length(SquidProfile(i).depthA_Thermocline10),1), SquidProfile(i).depthA_Thermocline10, '.k', 'MarkerSize', 15);
%     plot(repmat(i, length(SquidProfile(i).depthD_Thermocline10),1), SquidProfile(i).depthD_Thermocline10, '.m', 'MarkerSize', 15);
%     set(gca, 'YDir', 'reverse');
% end
% title([tagName ' 10° thermocline depths: day (black) and night (pink) '], 'FontSize', 16, 'FontWeight', 'bold')
% xlabel('Deployment Day ', 'FontSize', 14, 'FontWeight', 'bold')
% ylabel('Depth (m)','FontSize', 14, 'FontWeight', 'bold')
% set(h1, 'fontsize', 12, 'fontweight', 'bold')
% FigName = [tagName '_Thermocline10b.pdf'];
% annotate_JS(Mfilename, gcf, FigName);
% orient landscape % save
% print('-dpdf', [dir1 FigName]);

%% plot min and max temp/depth.

% subplot(1,2,1) % didn't fix this up on Mar 17 2011; not the most useful analysis.
% h(2) = area(DeployDates, Temptrix(:,2)); % max
% hold on
% h(1) = area(DeployDates, Temptrix(:,1)); % min
% set(h(2),'FaceColor',[0.76 0.87 0.78])
% set(h(1),'FaceColor',[1 1 1])
% set(gca,'Layer','top')
% plot(DeployDates, Temptrix(:,1:2), 'k','LineWidth',2)
% set(gca, 'fontsize', 12, 'fontweight', 'bold')
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
% set(gca, 'fontsize', 12, 'fontweight', 'bold')
% ylabel('Depth [m]', 'FontSize',20, 'FontWeight', 'bold')
% set(gca,'XTick', DeployDates)
% datetick('x','mm-dd','keepticks')
% [ax,h(1)] = suplabel('Date ','x');
% [ax,h(2)] = suplabel('','y');
% set(h,'FontSize',20, 'FontWeight', 'bold')
% key = {[tagName ' Daily Min and Max Temp and Depth, Series Data']; [datestr(min(DeployDates),22) ' - '  datestr(max(DeployDates),22)]};
% [h3] = suptitle(key);
% set(h3,'FontSize',24, 'FontWeight', 'bold')
% FigName = [tagNumS 'MinMaxTD.pdf'];
% annotate_JS(Mfilename, gcf, FigName)
% orient landscape
% print('-dpdf', [dir1 FigName]);
%
% % $fixfig mlc
%
% clear ax

%% datetime vs depth with temp as color

if pf
    figure(500); clf('reset');
    set(gcf,'Position',[12   213   842   688])
    scatter(Seriesdate,DepthTag*-1,20,TempTag,'filled');
    title([tagName ' Seriesdata2'])
    ylabel('Depth (m)')
    xlabel('Date (PST)')
    Hcbar = colorbar('location','SouthOutside', 'LineWidth', 1);
    caxis(colorsLoc);
    xlabel(Hcbar,'Temperature (\circC)');
    datetick('x',6)
    if exist('DayNight_boxes.m','file') %if the function is present the figure will have night boxes
        [x,y] = DayNight_boxes(Seriesdate,sr/24,ss/24,get(gca,'Ylim'),0);
        figure(500); clf('reset');
        set(gcf,'Position',[12   213   842   688])
        fill(x,y,[0.8 0.8 0.8], 'EdgeColor','none')
        hold on
        scatter(timeN,DepthN*-1,20,TempN,'filled');
        scatter(timeD,DepthD*-1,20,TempD,'filled');
        %     plot(Seriesdate(iNight),DepthTag(iNight)*-1, '.k', 'MarkerSize', 15); % to test for Day/Night accuracy
        %     plot(Seriesdate(iDay),DepthTag(iDay)*-1, '.r', 'MarkerSize', 15);
        hold off
        title([tagName ' ' buflab],'FontSize',20, 'FontWeight', 'bold')
        ylabel('Depth (m)', 'FontSize',16, 'FontWeight', 'bold')
        xlabel('Date (PST)', 'FontSize',16, 'FontWeight', 'bold')
        set(gca, 'XLim', [min(Seriesdate) max(Seriesdate)])
        set(gca, 'fontsize', 12, 'fontweight', 'bold')
        Hcbar = colorbar('location','SouthOutside', 'LineWidth', 1);
        caxis(colorsLoc); % defined in importTagData_csv.m
        xlabel(Hcbar,'Temperature (\circC)', 'FontSize',16, 'FontWeight', 'bold');
        set(gca,'XTick', DeployDates+ datenum(0, 0, 0, 8, 0, 0))
        datetick('x','mm/dd', 'keeplimits','keepticks')
    end
    FigName = [tagName '_TimeSeries_TimeDepth_Temp' buflab '.pdf'];
    annotate_JS(Mfilename, gcf, FigName);
%         orient landscape % save
%         print('-dpdf', [dir1 FigName]);
    
% and for black and white:
plotSeriesDataPublication

    
    %% vertical migrations, baby!
    
    Mfilename = mfilename;
    
    SeriesdateNtimes = datenum(SeriesdateVtimes);
    DepthDescenti = DepthDescent(iDescent);
    TempDescenti = TempDescent(iDescent);
    SeriesdateNtimesiD = SeriesdateNtimes(iDescent);
    DepthAscenti = DepthAscent(iAscent);
    TempAscenti = TempAscent(iAscent);
    SeriesdateNtimesiA = SeriesdateNtimes(iAscent);
    
    [x,y] = DayNight_boxes(Seriesdate,sr/24,ss/24,get(gca,'Ylim'),0);
    figure
    set(gcf,'Position',[12   213   842   688])
    fill(x,y,[0.8 0.8 0.8], 'EdgeColor','none')
    hold on
    plot(timeN,DepthN*-1,'.k', timeD,DepthD*-1, '.k');
    plot(timeAscent, DepthAscent*-1, '.m', timeDescent, DepthDescent*-1, '.g')
    hold off
    title([tagName 'Vertical Migrations'],'FontSize',20, 'FontWeight', 'bold')
    ylabel('Depth (m)  ', 'FontSize',16, 'FontWeight', 'bold')
    xlabel('Date (PST)  ', 'FontSize',16, 'FontWeight', 'bold')
    set(gca, 'XLim', [min(Seriesdate) max(Seriesdate)])
    set(gca, 'fontsize', 12, 'fontweight', 'bold')
    set(gca,'XTick', DeployDates+ datenum(0, 0, 0, 8, 0, 0))
    datetick('x','mm/dd', 'keeplimits','keepticks')
    FigName = [tagName '_TimeSeries_TimeDepth_TempVertMigrations.pdf'];
    annotate_JS(Mfilename, gcf, FigName);
    %     orient landscape % save
    %     print('-dpdf', [dir1 FigName]);
    
    % time vs. depth for each descent
    figure;
    set(gcf,'Position',[12   213   842   688])
    plot(SeriesdateNtimesiD-datenum(0,0,0,8,0,0),DepthDescenti, '.', 'Color', [0 0.5020 0]);
    ylabel('Depth [m]', 'FontSize',20, 'FontWeight', 'bold')
    xlabel('Time, PST', 'FontSize',20, 'FontWeight', 'bold')
    title([tagName ' Depth v. Time, Descent. SeriesData'], 'FontSize',20, 'FontWeight', 'bold')
    datetick('x', 13);
    set(gca, 'YDir', 'reverse');
    grid on
    FigName = [tagName '_Profile_DepthTimeD.pdf'];
    annotate_JS(Mfilename, gcf, FigName)
    %     orient landscape % save
    %     print('-dpdf', FigName);
    
    % time vs. depth for each Ascent
    figure;
    set(gcf,'Position',[12   213   842   688])
    plot(SeriesdateNtimesiA-datenum(0,0,0,8,0,0),DepthAscenti, '.', 'Color', [0 0.5020 0]);
    ylabel('Depth [m]', 'FontSize',20, 'FontWeight', 'bold')
    xlabel('Time, PST', 'FontSize',20, 'FontWeight', 'bold')
    title([tagName ' Depth v. Time, Ascent. SeriesData'], 'FontSize',20, 'FontWeight', 'bold')
    datetick('x', 13);
    set(gca, 'YDir', 'reverse');
    FigName = [tagName '_Profile_DepthTimeA.pdf'];
    annotate_JS(Mfilename, gcf, FigName)
    %     orient landscape % save
    %     print('-dpdf', FigName);
    
    
    % now try to figure out rates?
    % % time vs. depth derivative  for each descent
    
    figure; hold on
    for i = 1:length(DeployDays)-1
        dlog = DEPLOYINDEX == DeployDays(i);
        diD = dlog + iDescent;
        diDlog = diD == 2;
        depthder = diff(DepthTag(diDlog));
        depthder = [depthder; NaN];
        plot(SeriesdateNtimes(diDlog)-datenum(0,0,0,8,0,0), depthder, 'ok', 'MarkerSize', 5) % Figure out how to cycle color
        if i >= 3
            plot(SeriesdateNtimes(diDlog)-datenum(0,0,0,8,0,0), depthder, 'og', 'MarkerSize', 5)
        end
    end
    ylabel({'Depth Change (diff(Depth), in m/s)'; 'green are deep dives on 12/8, 12/9'} , 'FontSize',20, 'FontWeight', 'bold')
    xlabel('Time, PST', 'FontSize',20, 'FontWeight', 'bold')
    title([tagName ' Depth Derivative v. Time, Ascent. SeriesData'], 'FontSize',20, 'FontWeight', 'bold')
    datetick('x', 13);
    FigName = [tagName '_Profile_DepthDerivTimeA.pdf'];
    annotate_JS(Mfilename, gcf, FigName)
    %     orient landscape % save
    %     print('-dpdf', FigName);
    
end

%% Deep Dives

deepDiveVelo = input('Would you like deep dive analysis? (1 = Y, 0 = N): ');
if deepDiveVelo % for deep dives only. Defined in manageTagData.m. Code taken from calcOMZtimes.m
    
    %     LB = 1200; % meters
    %
    if TagLocation == 1
        D_UB = 484.5; % median OMZ depth
%         if tagID == num2str(1)
%             DepthTag = DepthTag*-1;
%         end
        
    else
        switch tagID
            case num2str(4)
                D_UB = 204.6;% median OMZ depths
            case num2str(5)
                D_UB = 204.6;
            case num2str(6)
                D_UB = 247.8;
            case num2str(7)
                D_UB = 370.6;
        end
        %         DepthTag = DepthTag*-1;
    end
    
    
    close all
    xyData = []; % x,y input for each tag, make it long enough so can accommodate all events
    figure; hold on
    plot(Seriesdate,abs(DepthTag));
    plot(Seriesdate,repmat(D_UB,length(Seriesdate),1), '--r')
    set(gca, 'YDir', 'reverse')
    grid minor
    grid on
    title({'ginput for deep dives: '; ...
        ['identify (1) start of dive, (2) OMZ upper boundary at ' num2str(D_UB) ' m, ']; ...
        '(3) OMZ lower boundary at 1200 m, (4) max, (5) OMZ LB, (6) OMZ UB, (7) end of dive'})
    [t,d] = ginput;
    d = abs(d);
    close all
    
    deltaT = zeros(length(t)-1,1); % in minutes
    deltaD = zeros(length(t)-1,1);
%     ind1 = find(Seriesdate == t(1));
    
    for k = 1:length(t)-1
        tt = t(k+1)-t(k); 
        dd = d(k)-d(k+1);
        deltaT(k) = tt;
        deltaD(k) = dd;
    end
    
    deltaTa = datevec(deltaT);
    deltaT_H = deltaTa(:,4) + deltaTa(:,5)/60 + deltaTa(:,6)/3600; % this is just fact checking with timing

%     bb = deltaD./deltaT_H;
%     bb/3600
    
%     for pp = 1:length(deltaT)
%         ppp = deltaD(pp)/deltaT_H(pp)
%         ppp/3600
%     end
    
    % to put in .xls sheet TaggingLog.xls
    d;
    t;
    t_string = datestr(t, 'HH:MM:SS');
    deltaD;
    deltaT;
    deltaT_hr = deltaT*24;
    veloDT_m_hr = deltaD./deltaT_hr;
    veloDT_m_min = deltaD./(deltaT*24*60);
    veloDT_m_sec = deltaD./(deltaT*24*60*60);
    
end

% a = d(2)-d(1)
% 
% t(2)-t(1)
% 
% MIN = a/19 % minutes
% MIN/60
% 
% deltaT_H./deltaD
% 
% aa=d(4)-d(3)
% t(4)-t(3)
% MIN2 = a/38
% MIN2/60

% save info on deep dives. 

% save([ tagName '_DeepDiveInfo.mat'] , 'd', 't', 't_string', 'veloDT_m_hr', 'veloDT_m_sec') 

%%

% % rate?
% figure; hold on
% for i = 1:length(DeployDays)
%     dlogA = DEPLOYINDEX == DeployDays(i);
%     diA = dlogA + iAscent;
%     diAlog = diA == 2;
%
%     plot(1:sum(diAlog), DepthTag(diAlog)./SeriesdateNtimes(diAlog), '.') % Figure out how to cycle color
% end
% title([datestr(DeployDates(i),6) ', ' num2str(ceil(DeployHours(i))) ' hours'], 'FontSize',12, 'FontWeight', 'bold')
% xlabel('Rate m/75sec, is that right?', 'FontSize',20, 'FontWeight', 'bold')
% set(gca, 'YDir', 'reverse');

%% save all figures Concatenated

% Mfilename = mfilename; % just to reinstate
% Fname = ['_' tagName '_CombinedFigs_' Mfilename];
% FigHs = unique(findobj(0, '-depth',1, 'type','figure')); %get all open figure handles
% for F = 1:length(FigHs)
%     if F == 1
%         print ( '-dpsc2', [Fname '.ps'], ['-f' num2str(FigHs(F))])%,'-painters']) %create postscript file of figure
%     else
%         print ( '-dpsc2', [Fname '.ps'], '-append', ['-f' num2str(FigHs(F))])%,'-painters']) %append next figure to postscript file
%     end
% end
% ps2pdf('psfile', [Fname '.ps'], 'pdffile',  [Fname '.pdf'], 'deletepsfile', 1) %ps2pdf.m MATLAB Central, convert postscript file to pdf

% NOTE: to just save little bits, but not everything, when you're making
% individual figures, put the figure number into the vector FigHs.

%%

disp('Completed plotSeriesData.m')
% ===== EOF [plotSeriesData.m] ======
