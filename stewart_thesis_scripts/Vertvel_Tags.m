% Vertvel_Tags.m
%
% Outputs:
%
% Outside Functions Called:
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 04-Mar-2011 12:34:01
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : Vertvel_Tags.m

Mfilename = mfilename;

% Change DepthRaw to DepthTag if you want to run 83051 with 1 sec
% sampling...

%% Calculate Horizontal distance

if tagNumD == 83051 % this is done for all tags in manageTagData_CTD.m 
    [dist_deg,az] = distance(36.5822, -122.0569333, 32.09194444, -118.4463889); % 83051 deploy and popup
    dist_km = deg2km(dist_deg);
    DeployHoursTot = sum(DeployHours);
    dist_km_hr = dist_km/DeployHoursTot;
    dist_km_da = dist_km_hr*24;
end


%%
if TagRecovered
    Rate = 1; % archival data at 1 second
    VertbinWidth = 2;
    VertbinCtrs = 0:VertbinWidth:100;
    YmaxCount = 40000;
    DEPTH = DepthTag;
    TIMEV = SeriesdateV;
else
    Rate = 75; % uploaded timeseries data at 75 seconds
    VertbinWidth = 2;
    VertbinCtrs = 0:VertbinWidth:60;
    YmaxCount = 800;
    DEPTH = DepthRaw;
    TIMEV = DateRawV;
end

%%

Depth_sm = sgolayfilt(DEPTH,6,33); %smooths trace while keeping large peaks % Depth raw makes all tags at 75 sec
Vertvel = diff(Depth_sm)./Rate;
Vertvel_m_sec =  [diff(Depth_sm)./Rate; NaN];%1st derivative
Vertvel_m_sec(Vertvel_m_sec < 0) = Vertvel_m_sec(Vertvel_m_sec < 0)*-1; % get rid of negative velocities, just looking for absolute value
summary(Vertvel_m_sec)

m_sec_to_km_da = 60*60*24/1000;
km_day_to__m_sec = 1000/(60*60*24);
Vertvel_km_day = Vertvel_m_sec*m_sec_to_km_da;
summary(Vertvel_km_day)

%% plot counts

% figure out vertvel for each day
IQRquant = [.05 .5 .95];
Vertvel_km_day_Max_Quant = zeros(length(DeployDates),4);

Count = 0;
for i = 1:length(DeployDates)
    dlog = TIMEV(:,3) == DeployDatesV(i,3);
    Vertvel_km_day_i = Vertvel_km_day(dlog);
    Vertvel_km_day_ibinned = hist(Vertvel_km_day_i,VertbinCtrs); %count within each bin
    Vertvel_km_day_iPerc = ((Vertvel_km_day_ibinned/sum(Vertvel_km_day_ibinned))*100)';
    %     Vertvel_m_secSummary = summary(Vertvel_km_day_i);
    Vertvel_km_day_Max_Quant(i,:) = [quantile(Vertvel_km_day_i, IQRquant) max(Vertvel_km_day_i)];
    
    Count = Count + 1;
    fignum = floor((Count - 1)./numberofpanels)+1; % set the number for 6-panel figures
    fignum2 = fignum+90;
    panelnum = rem(Count - 1, numberofpanels)+1; % rem takes the remainder after division
    if panelnum == 1
        figure(fignum2); clf;
    end
    subplot(ceil(numberofpanels/numcols), numcols, panelnum); hold on;
    hist(Vertvel_km_day_i, VertbinCtrs);
    counts = hist(Vertvel_km_day_iPerc, VertbinCtrs);
    set(gca,'XLim', [min(VertbinCtrs) max(VertbinCtrs)], 'YLim', [0 YmaxCount]);
    text(max(VertbinCtrs)*.3, max(YmaxCount)*.7, {['Vel max = ' num2str(ceil(max(Vertvel_km_day_i))) ' '];...
        ['Vel at ' num2str(IQRquant*100) '%: ']; ['       ' num2str(quantile(Vertvel_km_day_i, IQRquant)) ' ']}, 'fontsize', 12, 'fontweight', 'bold')
    title([datestr(DEPLOYINDEXuniN(i),6) ', ' num2str(ceil(DeployHours(i))) ' hours'], 'FontSize',12, 'FontWeight', 'bold') % error if spans 2 months--fix if need be using DeployMonth
    
    if panelnum == 6
        key = {[tagName ' Daily Vertical Velocities ' num2str(fignum)]};
        [h3] = suptitle(key);
        [ax,h(1)] = suplabel('Velocity (km/day) ','x');
        [ax,h(2)] = suplabel('Count ','y');
        set(h3,'FontSize',24, 'FontWeight', 'bold')
        set(h,'FontSize',20, 'FontWeight', 'bold')
        FigName = [tagNumS 'Series_daily_Vertvel' num2str(fignum) '.pdf'];
        annotate_JS(Mfilename, gcf, FigName)
%         orient landscape % save
%         print('-dpdf', [dir1 FigName]);
    end
end
if panelnum < 6 % catch for when there aren't 6 panels on a figure
    key = {[tagName ' Daily Vertical Velocities ' num2str(fignum)]};
    [h3] = suptitle(key);
    [ax,h(1)] = suplabel('Velocity (km/day) ','x');
    [ax,h(2)] = suplabel('Count ','y');
    set(h3,'FontSize',24, 'FontWeight', 'bold')
    set(h,'FontSize',20, 'FontWeight', 'bold')
    FigName = [tagNumS 'Series_daily_Vertvel' num2str(fignum) '.pdf'];
    annotate_JS(Mfilename, gcf, FigName)
%     orient landscape % save
%     print('-dpdf', [dir1 FigName]);
end

Vertvel_km_day_Max_Quant

%% %% plot percents. IN PROGRESS

% % figure out vertvel for each day
% VertbinWidth = 2;
% VertbinCtrs = 0:VertbinWidth:60;
% IQRperc = 90;
% Count = 0;
% Vertvel_km_day_i_Tag = zeros(length(DeployDates),1);
% for i = 1:length(DeployDates)
%     dlog = SeriesdateV(:,3) == DeployDatesV(i,3);
%     Vertvel_km_day_i = Vertvel_km_day(dlog);
%     Vertvel_km_day_ibinned = hist(Vertvel_km_day_i,VertbinCtrs); %count within each bin
%     Vertvel_km_day_iPerc = ((Vertvel_km_day_ibinned/sum(Vertvel_km_day_ibinned))*100)';
% %     Vertvel_m_secSummary = summary(Vertvel_km_day_i);
% %     Vertvel_km_day_i_Tag(i) = Vertvel_km_day_i;
% 
%     Count = Count + 1;
%     fignum = floor((Count - 1)./numberofpanels)+1; % set the number for 6-panel figures
%     fignum2 = fignum+70;
%     panelnum = rem(Count - 1, numberofpanels)+1; % rem takes the remainder after division
%     if panelnum == 1
%         figure(fignum2); clf;
%     end
%     subplot(ceil(numberofpanels/numcols), numcols, panelnum); hold on;
%     hist(Vertvel_km_day_iPerc);%, VertbinCtrs);
%     counts = hist(Vertvel_km_day_iPerc, VertbinCtrs);
% 
% %
% %     h = get(gca,'child');
% % set(h,'FaceColor',[.98 .98 .98],'EdgeColor',[.94 .94 .94]);
% % % counts = hist(Vertvel_km_day_i,VertbinCtrs);
% % hold on
% % % plot(VertbinCtrs,counts,'o');
% % % hold off
% %
% %     Vertvel_km_day_i(isnan(Vertvel_km_day_i)) = [];
% %     paramEsts = wblfit(Vertvel_km_day_i);
% %     n = length(Vertvel_km_day_i);
% %     prob = counts / (n * VertbinWidth);
% %     bar(VertbinCtrs,prob,'hist');
% %     h = get(gca,'child');
% % %     set(h,'FaceColor',[.9 .9 .9]);
% % %     xlabel('Time to Failure'); ylabel('Probability Density'); );
% %     xgrid = linspace(0,20,100);
% %     pdfEst = wblpdf(xgrid,paramEsts(1),paramEsts(2));
% %     hold on
% %     plot(xgrid,pdfEst)
% 
% 
% 
% %     set(gca,'XLim', [min(VertbinCtrs) max(VertbinCtrs)], 'YLim', [0 1]);
%     text(max(VertbinCtrs)*.3, max(YmaxCount)*.7, {['Vel max = ' num2str(ceil(max(Vertvel_km_day_i))) ' '];...
%         ['Vel ' num2str(IQRperc) '% IQR = ']; ['       ' num2str(quantile(Vertvel_km_day_i, [.05 .95])) ' ']}, 'fontsize', 12, 'fontweight', 'bold')
%     title([datestr(datelabeltrix(i),6) ', ' num2str(ceil(DeployHours(i))) ' hours'], 'FontSize',12, 'FontWeight', 'bold') % error if spans 2 months--fix if need be using DeployMonth
% 
%     if panelnum == 6
%         key = {[tagName ' Daily Vertical Velocities ' num2str(fignum)]};
%         [h3] = suptitle(key);
%         [ax,h(1)] = suplabel('Velocity (km/day) ','x');
%         [ax,h(2)] = suplabel('Percent ','y');
%         set(h3,'FontSize',24, 'FontWeight', 'bold')
%         set(h,'FontSize',20, 'FontWeight', 'bold')
%         FigName = [tagNumS 'Series_daily_Vertvel' num2str(fignum) '.pdf'];
%         annotate_JS(Mfilename, gcf, FigName)
%         orient landscape % save
%         print('-dpdf', [dir1 FigName]);
%     end
% end
% if panelnum < 6 % catch for when there aren't 6 panels on a figure
%     key = {[tagName ' Daily Vertical Velocities ' num2str(fignum)]};
%     [h3] = suptitle(key);
%     [ax,h(1)] = suplabel('Velocity (km/day) ','x');
%     [ax,h(2)] = suplabel('Percent ','y');
%     set(h3,'FontSize',24, 'FontWeight', 'bold')
%     set(h,'FontSize',20, 'FontWeight', 'bold')
%     FigName = [tagNumS 'Series_daily_Vertvel' num2str(fignum) '.pdf'];
%     annotate_JS(Mfilename, gcf, FigName)
% %     orient landscape % save
% %     print('-dpdf', [dir1 FigName]);
% end
% 

%%

disp('Completed Vertvel_Tags.m')
% ===== EOF [Vertvel_Tags.m] ======
