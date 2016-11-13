% tagDailyHistos_ShallowPDT.m
% called from processPDTData.m % J. Stewart
%
% this script will run and save daily histo data for each tag, to be used
% with SST Hinke analysis.Only uses PDT data.
%
%
% called from processSeriesData.m by J. Stewart
%
% Outside Functions Called:
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 25-Feb-2011 11:02:37
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : tagDailyHistos_ShallowPDT.m

Mfilename = mfilename;

ShallowCutoff = 50%5:5:50;
for sc = 1:length(ShallowCutoff)
    
    %% setup lifted from tagDailyHistos_Shallow.m
    
    DepthShallowlog = DepthTag <= ShallowCutoff(sc); % ShallowCutoff(sc) defined in manageTagData
    DepthShallow = DepthTag(DepthShallowlog,1);
    TempShallow = TempTag(DepthShallowlog,1); %TempShallow = TempTagMin(DepthShallowlog,1); % actually, there's little difference between TempTagMin and TempTag.
    DatesTagShallow = DEPLOYINDEXN(DepthShallowlog,1);
    DatesTagShallowV = datevec(DatesTagShallow);
    DepthTempTag = [DepthTag, TempTag];
    
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
    
    %% temp vs. depth daily subplots
    Count = 0;
    numberofpanels = 6;
    numcols = 3;
    
    maxtemp = 16;
    mintemp = 9;
    maxcount = 30;
    binset = 0.1;
    binincr = mintemp:binset:maxtemp; % because hist works with centers
    DepthTempTagDay = [];
    Count = 0;
    HistCounts = []; % summary is not the way to go because it's not normal at all.
    ShallowDailyTemp = [];
    for i = 1:length(DeployDates)
        ilog = DatesTagShallow == DeployDates(i);
        tlog = DEPLOYINDEXN == DeployDates(i);
        
        TagLength = length(TempTag) % total amount of data counts for all tags
        TagShallowLength = length(TempShallow(ilog))
        TagShallowPerc = (TagShallowLength/TagLength)*100
        
        DepthTempTagDayA = DepthTempTag(tlog,:);
        DepthTempTagDayB = sortrows(DepthTempTagDayA,1); % sort by first row (ie depth)
        DepthTempTagDayC = DepthTempTagDayB(~isnan(DepthTempTagDayB(:,2)),:); % troubleshoot for days where all nans
        if i ~= length(DeployDates)
            DepthTempTagDay = [DepthTempTagDay DepthTempTagDayC(1:5,:)];
            %         else
            %             DepthTempTagDay = [DepthTempTagDay DepthTempTagDayC]; % last day has only 9 hits
        end
        Count = Count + 1;
        fignum = floor((Count - 1)./numberofpanels)+1; % set the number for 6-panel figures
        fignum2 = fignum+100;
        panelnum = rem(Count - 1, numberofpanels)+1; % rem takes the remainder after division
        if panelnum == 1
            figure(fignum2); clf;
        end
        subplot(ceil(numberofpanels/numcols), numcols, panelnum); hold on;
        histtemp = hist(TempShallow(ilog), binincr);
        hist(TempShallow(ilog), binincr)
        h = findobj(gca, 'Type', 'patch');
        set(h, 'FaceColor', [0.701960802078247 0.780392169952393 1]);
        if size(histtemp,2) < size(histtemp,1) % if it's a row not column vector
            histtemp = histtemp';
        end
        HistCounts = [HistCounts; histtemp];
        ShallowDailyTemp = [ShallowDailyTemp; repmat(DeployDays(i),sum(ilog),1) TempShallow(ilog)];
        %         for j = 1:length(binincr)
        %             text(binincr(j), maxcount-2, num2str(histtemp(j)))
        %         end
        set(gca, 'XLim', [mintemp maxtemp], 'YLim', [0 maxcount])
        title([datestr(DeployDates(i),6)], 'FontSize',12, 'FontWeight', 'bold')
        if panelnum == 6
            key = {[tagName ' Daily Temp Shallow at ' num2str(ShallowCutoff(sc)) 'm: '];...
                [num2str(ceil(TagShallowPerc)) '% of tag data, ' num2str(binset) '\circC bins. Figure ' num2str(fignum)]};
            [h3] = suptitle(key);
            [ax,h(1)] = suplabel('Depth [m] ','x'); %('Day of the Year ','x');
            [ax,h(2)] = suplabel('Temperature [\circC] ','y');
            set(h3,'FontSize',24, 'FontWeight', 'bold')
            set(h,'FontSize',20, 'FontWeight', 'bold')
            FigName = [tagName '_PDT_tempbinsDailyShallow' num2str(ShallowCutoff(sc)) 'm' num2str(fignum) '.pdf'];
            annotate_JS(Mfilename, gcf, FigName)
%             orient landscape % save
%             print('-dpdf', [dir1 FigName]);
        end
    end
    if panelnum < 6 % catch for when there aren't 6 panels on a figure
        key = {[tagName ' Daily Temp Shallow at ' num2str(ShallowCutoff(sc)) 'm: '];...
            [num2str(ceil(TagShallowPerc)) '% of tag data, ' num2str(binset) '\circC bins. Figure ' num2str(fignum)]};
        [h3] = suptitle(key);
        [ax,h(1)] = suplabel('Temperature [\circC] ','x');
        [ax,h(2)] = suplabel('Depth [m] ','y');
        set(h3,'FontSize',24, 'FontWeight', 'bold')
        set(h,'FontSize',20, 'FontWeight', 'bold')
        FigName = [tagName '_PDT_tempbinsDailyShallow' num2str(ShallowCutoff(sc)) 'm' num2str(fignum) '.pdf'];
        annotate_JS(Mfilename, gcf, FigName)
%         orient landscape % save
%         print('-dpdf', [dir1 FigName]);
    end
    
end
%% make a table of the hiscounts for each day

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
% % orient landscape % save
% % print('-dpdf', [dir1 FigName]);

%% Make matrices for each tag based on the tables, and save them as .mat files

DailyShallowHistCounts = HistCounts;

% now save it:
% save -mat DailyShallow_83046_45.mat DailyShallowHistCounts binincr ShallowDailyTemp DepthTagMax % make sure it's labeled for the proper depth and tag num
save -mat DailyShallowest5_83046.mat DepthTempTagDay DailyShallowHistCounts binincr ShallowDailyTemp DepthTagMax % make sure it's labeled for the proper depth and tag num

%%

disp('Completed tagDailyHistos_ShallowPDT.m')
% ===== EOF [tagDailyHistos_ShallowPDT.m] ======
