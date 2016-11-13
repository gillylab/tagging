% plotVertMigration.m
%  
% Called From: processSeriesData.m % J. Stewart
  
% Description: 
%  
% Outside Functions Called: 
% findDielMigrationJS.m % J. Stewart
%  
% AUTHOR    : J. Stewart jules32@gmail.com 
% DATE      : 02-Jun-2011 09:52:49  
% Revision  : 1.00  
% DEVELOPED : 7.9.0.529 (R2009b) OSX 
% FILENAME  : plotVertMigration.m 

Mfilename = mfilename;  
  
%% plot this to identify the hour cutoffs for each tag (bufferA, bufferD in processSeriesData).

[LegLabel] = makeLabels(1:length(DeployDays), ' day'); % A. Booth

figure; hold on
for i = 1:length(DeployDays) % have to plot twice so legend shows up
    plot(DescentDatetrix(:,i), DescentDepthtrix(:,i), '.', 'Color', Col(i,:)) % Figure out how to cycle color
end
legend(LegLabel,'Location', 'SouthWest');
set(gca, 'YDir', 'reverse');
[x,y] = DayNight_boxes(Seriesdate,sr/24,ss/24,get(gca,'Ylim'),0);
set(gcf,'Position',[12   213   842   688])
fill(x,y,[0.8 0.8 0.8], 'EdgeColor','none')
hold on
plot(timeN,DepthN,'.k',timeD,DepthD,'.k')
for i = 1:length(DeployDays)
    plot(DescentDatetrix(:,i), DescentDepthtrix(:,i), '.m')%, 'Color', Col(i,:)) % Figure out how to cycle color
    plot(AscentDatetrix(:,i), AscentDepthtrix(:,i), '.c')%,  'Color', Col(i,:)) % Figure out
end
title([tagName ' Diel Migration Markings   '], 'FontSize',18, 'FontWeight', 'bold')
xlabel('Time ', 'FontSize',14, 'FontWeight', 'bold')
ylabel('Depth (m)', 'FontSize',14, 'FontWeight', 'bold')

set(gca,'XTick', DeployDates+datenum(0, 0, 0, 8, 0, 0))
datetick('x','mm/dd', 'keeplimits','keepticks')
FigName = [tagName '_SeriesDailyMigrations.pdf'];
annotate_JS(Mfilename, gcf, FigName)
% orient landscape % save
% print('-dpdf', [dir1 FigName]);

%% Plot Daily Descents all ontop of each other
ifplot = 1;
if ifplot
    DescentDatetrixVtimes = [];
    DeployS = num2str(1:length(DeployDays));
    
    for i = 1:length(DeployDays)
        Dtemp1 = datevec(DescentDatetrix(:,i));
        Dtemp1a = datenum(0,0,0,Dtemp1(:,4),Dtemp1(:,5),Dtemp1(:,6));
        DescentDatetrixVtimes = [DescentDatetrixVtimes Dtemp1a];
    end
    
    minDD = min(min(DescentDatetrixVtimes));%+datenum(0, 0, 0, 8, 0, 0);
    maxDD = max(max(DescentDatetrixVtimes));%+datenum(0, 0, 0, 8, 0, 0);
    XlabA = minDD:0.05:maxDD;
    Xlab = XlabA+datenum(0, 0, 0, 8, 0, 0);
    
    figure; hold on
    for i = 1:length(DeployDays) % have to plot this twice so the legend works
        plot(DescentDatetrixVtimes(:,i), DescentDepthtrix(:,i), '.', 'Color', Col(i,:)) %-datenum(0,0,0,8,0,0)
    end
    legend(LegLabel,'Location', 'SouthWest');
    [x,y] = DayNight_boxes(DescentDatetrixVtimes(:,3),sr/24,ss/24,get(gca,'Ylim'),0);
    set(gcf,'Position',[12   213   842   688])
    fill(x,y,[0.8 0.8 0.8], 'EdgeColor','none')
    for i = 1:length(DeployDays)
        plot(DescentDatetrixVtimes(:,i), DescentDepthtrix(:,i), 'Color', Col(i,:)) %-datenum(0,0,0,8,0,0)
    end
    title([tagName ' Diel Migration Descents   '], 'FontSize',18, 'FontWeight', 'bold')
    xlabel('Time, GMT (shading is PST)', 'FontSize',14, 'FontWeight', 'bold')
    ylabel('Depth (m)', 'FontSize',14, 'FontWeight', 'bold')
    set(gca, 'YDir', 'reverse');
    set(gca,'XLim', [minDD maxDD])
    set(gca,'XTickLabel', num2str(Xlab))
    datetick('x',15,'keepticks')
    FigName = [tagName '_SeriesDailyDescents.pdf'];
    annotate_JS(Mfilename, gcf, FigName)
%     orient landscape % save
%     print('-dpdf', [dir1 FigName]);
    
    %% Plot Daily Ascents all ontop of each other
    
    AscentDatetrixVtimes = [];
    
    for i = 1:length(DeployDays)
        Atemp1 = datevec(AscentDatetrix(:,i));
        Atemp1a = datenum(0,0,0,Atemp1(:,4),Atemp1(:,5),Atemp1(:,6));
        AscentDatetrixVtimes = [AscentDatetrixVtimes Atemp1a];
    end
    
    minAD = min(min(AscentDatetrixVtimes));%+datenum(0, 0, 0, 8, 0, 0);
    maxAD = max(max(AscentDatetrixVtimes));%+datenum(0, 0, 0, 8, 0, 0);
    XlabA = minAD:0.05:maxAD;
    Xlab = XlabA+datenum(0, 0, 0, 8, 0, 0);
    
    % plot
    figure; hold on
    for i = 1:length(DeployDays) % have to plot this twice so the legend works
        plot(AscentDatetrixVtimes(:,i), AscentDepthtrix(:,i), '.', 'Color', Col(i,:)) %-datenum(0,0,0,8,0,0)
    end
    legend(LegLabel,'Location', 'SouthEast');
    [x,y] = DayNight_boxes(AscentDatetrixVtimes(:,3),sr/24,ss/24,get(gca,'Ylim'),0);
    set(gcf,'Position',[12   213   842   688])
    fill(x,y,[0.8 0.8 0.8], 'EdgeColor','none')
    for i = 1:length(DeployDays)
        plot(AscentDatetrixVtimes(:,i), AscentDepthtrix(:,i), 'Color', Col(i,:)) %-datenum(0,0,0,8,0,0)
    end
    title([tagName ' Diel Migration Ascents   '], 'FontSize',18, 'FontWeight', 'bold')
    xlabel('Time, GMT (shading is PST)', 'FontSize',14, 'FontWeight', 'bold')
    ylabel('Depth (m)', 'FontSize',14, 'FontWeight', 'bold')
    set(gca, 'YDir', 'reverse');
    set(gca,'XLim', [minAD maxAD])
    set(gca,'XTickLabel', num2str(Xlab))
    datetick('x',15,'keepticks')
    FigName = [tagName '_SeriesDailyAscents.pdf'];
    annotate_JS(Mfilename, gcf, FigName)
%     orient landscape % save
%     print('-dpdf', [dir1 FigName]);
    
end

%% Ashley's approach, right now not working. 

% [DepthDA,DepthDD] = findDielMigrationJS(DepthTag,Seriesdate,max(DepthTag),min(DepthTag),ss,sr);

%% 
  
disp('Completed plotVertMigration.m') 
% ===== EOF [plotVertMigration.m] ======  
