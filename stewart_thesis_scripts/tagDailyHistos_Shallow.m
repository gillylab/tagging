% tagDailyHistos_Shallow.m
%
% this script will run and save daily shallow histo data for each tag, to be used with SST Hinke analysis.
% will save as .mat files. 
%
% called from processSeriesData.m by J. Stewart
%
% Outside Functions Called:
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 25-Feb-2011 11:02:37
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : tagDailyHistos_Shallow.m

Mfilename = mfilename;

%%
ShallowCutoff = 50; %10:10:50;
for sc = 1:length(ShallowCutoff)
    
    %% setup
    
    DepthShallowlog = DepthTag <= ShallowCutoff(sc); % ShallowCutoff(sc) defined in manageTagData
    DepthShallow = DepthTag(DepthShallowlog,1);
    TempShallow = TempTag(DepthShallowlog,1);
    DepthTempTag = [DepthTag, TempTag];
    
    DatesTagShallow = DEPLOYINDEXN(DepthShallowlog,1);
    DatesTagShallowV = datevec(DatesTagShallow);
    
    % percent of data
    TagLength = length(TempTag) % total amount of data counts for all tags
    TagShallowLength = length(TempShallow)
    TagShallowPerc = (TagShallowLength/TagLength)*100
    
    % set X limits
    shallowcutoff = 50; % can change this/get rid of after you decide on a ShallowCutoff(sc).
    scolog = DepthTag <= shallowcutoff; % ShallowCutoff(sc) defined in manageTagData
    scotemp = TempTag(scolog,1);
    mintemp = floor(min(scotemp))-1; % switch this back to TempTagShallow after it's settled.
    maxtemp = ceil(max(scotemp))+1;
    maxcount = 300; % figure out when to change that
    
    %% daily histos
    
    % temp
    binset = 0.1;
    binincr = mintemp:binset:maxtemp; % because hist works with centers
    Count = 0;
    HistCounts = []; 
    DepthTempTagDay = []; %50 shallowest depths w temps
    QuantileSet = [0.025, 0.05, 0.25 0.5, 0.75, 0.95, 0.975]; 
    TempShallowQuantile = []; 
    ShallowDailyTemp = [];
    for i = 1:length(DeployDates) 
        ilog = DatesTagShallow == DeployDates(i);
        tlog = DEPLOYINDEXN == DeployDates(i);
        DepthTempTagDayA = DepthTempTag(tlog,:);
        DepthTempTagDayB = sortrows(DepthTempTagDayA,1); % sort by first row (ie depth)
        DepthTempTagDayC = DepthTempTagDayB(~isnan(DepthTempTagDayB(:,2)),:); % troubleshoot for days where all nans
        DepthTempTagDay = [DepthTempTagDay DepthTempTagDayC(1:10,:)];
        
        Count = Count + 1;
        fignum = floor((Count - 1)./numberofpanels)+1; % set the number for 6-panel figures
        fignum2 = fignum+424;
        panelnum = rem(Count - 1, numberofpanels)+1; % rem takes the remainder after division
        if panelnum == 1
            figure(fignum2); clf;
        end
        subplot(ceil(numberofpanels/numcols), numcols, panelnum); hold on;
        hist(TempShallow(ilog), binincr);
        histtemp = hist(TempShallow(ilog), binincr);
        h = findobj(gca, 'Type', 'patch');
        set(h, 'FaceColor', [0.701960802078247 0.780392169952393 1]);
        if size(histtemp,1) > size(histtemp,2) % troubleshooting
            histtemp = histtemp';
        end
        HistCounts = [HistCounts; histtemp];
        tempquantile = quantile(TempShallow(ilog), QuantileSet);
        TempShallowQuantile = [TempShallowQuantile; tempquantile];
        ShallowDailyTemp = [ShallowDailyTemp; repmat(DeployDays(i),sum(ilog),1) TempShallow(ilog)];
        %     for j = 1:length(binincr)
        %         text(binincr(j), maxcount, num2str(histtemp(j)))
        %     end
        set(gca, 'XLim', [mintemp maxtemp], 'YLim', [0 maxcount])
        title([datestr(DeployDates(i),6) ', ' num2str(ceil(DeployHours(i))) ' hours'], 'FontSize',12, 'FontWeight', 'bold')
        if panelnum == 6
            key = {[tagName ' Daily Temp Shallow at ' num2str(ShallowCutoff(sc)) 'm: '];...
                [num2str(ceil(TagShallowPerc)) '% of tag data, ' num2str(binset) '\circC bins. Figure ' num2str(fignum)]}; 
            [h3] = suptitle(key);
            [ax,h(1)] = suplabel('Temperature (\circC bins) ','x'); %('Day of the Year ','x');
            [ax,h(2)] = suplabel('Count ','y');
            set(h3,'FontSize',24, 'FontWeight', 'bold')
            set(h,'FontSize',20, 'FontWeight', 'bold')
            FigName = [tagName '_Series_tempbinsDailyShallow' num2str(ShallowCutoff(sc)) 'm' num2str(fignum) '.pdf'];
            annotate_JS(Mfilename, gcf, FigName)
%             orient landscape % save
%             print('-dpdf', [dir1 FigName]);
        end
    end
    if panelnum < 6 % catch for when there aren't 6 panels on a figure
        key = {[tagName ' Daily Temp Shallow at ' num2str(ShallowCutoff(sc)) 'm: '];...
            [num2str(ceil(TagShallowPerc)) '% of tag data, ' num2str(binset) '\circC bins. Figure ' num2str(fignum)]};
        [h3] = suptitle(key);
        [ax,h(1)] = suplabel('Temperature (\circC bins) ','x');
        [ax,h(2)] = suplabel('Count ','y');
        set(h3,'FontSize',24, 'FontWeight', 'bold')
        set(h,'FontSize',20, 'FontWeight', 'bold')
        FigName = [tagName '_Series_tempbinsDailyShallow' num2str(ShallowCutoff(sc)) 'm' num2str(fignum) '.pdf'];
        annotate_JS(Mfilename, gcf, FigName)
%         orient landscape % save
%         print('-dpdf', [dir1 FigName]);
    end
    
    DailyShallowHistCounts = HistCounts;
    % OLD: filenameSave = ['DailyShallow_' tagNumS '_' num2str(ShallowCutoff(1)) '.mat']; % need to figure out how to make name for saving
%     filenameSave = ['DailyShallowest10_' tagNumS '.mat'];
%     eval (['save -mat ' filenameSave ' DepthTempTagDay DailyShallowHistCounts binincr ShallowDailyTemp DepthTagMax']) % DepthTagMax was made in processSeriesData.m
end

%% make a table of the hiscounts for each day. Runs, but not useful?

% ColHeader = [];
% for i = 1:length(binincr)
%     ColHeader = [ColHeader {binincr(i)}];
% end
% ColHeader = ['Temp Bin' ColHeader];
% RowHeader = datestr(DeployDates,6);
% 
% figure
% tbl = table(ColHeader, RowHeader, HistCounts);
% plot_table(tbl);
% [x,y] = getAxisInset( .1, .95 );
% pos = [360   502   670   420];
% set( gcf,'pos', pos );
% text(x,y, [tagName ' Temp Shallow: ' num2str(ShallowCutoff(sc)) ' m'], 'fontsize', 18,'fontweight', 'bold');
% FigName = [tagName '_TableTempShallow' num2str(ShallowCutoff(sc)) 'm.pdf'];
% annotate_JS(Mfilename, gcf, FigName)
% orient landscape % save
% print('-dpdf', [dir1 FigName]);

%%

disp('Completed tagDailyHistos_Shallow.m')
% ===== EOF [tagDailyHistos_Shallow.m] ======
