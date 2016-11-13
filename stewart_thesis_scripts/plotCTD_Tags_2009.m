% plotCTD_Tags_2009.m
%
% Called from manageTagData_CTD
%

% Section Titles: 
%         Determine which CTDs to use, based on date (user input)
%         Plot Oxygen data for chosen CTD, option to do this ontop of tag histogram
%         compare oxygens across CTD casts (also map of CTD locations)
%         compare temperatures across CTD casts
%         compare temperatures to oxygen across CTD casts
% Outside Functions Called:
% mapSqArea by J. Stewart: maps Monterey area of 2009 tags
% annotate_JS by J. Stewart
% b_ROVCTD_stats % by J. Stewart
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 01-Feb-2011 18:22:36
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : plotCTD_Tags_2009.m

Mfilename = mfilename;

%% Determine which CTDs to use, based on date (user input)

divenum_CTDtags2009 = [3403.1; 3417.1; 91.3; 3453.1; 101.3; 3481.1; 3514.1]; % after running this for a few months around squid tagging time, create this.
colorCTD = lines(length(divenum_CTDtags2009));

%% compare oxygens across CTD casts (also map of CTD locations)

figure
lavector = [];
lovector = [];
datvector = [];
oxyvector = [];
depvector = [];
for i = 1:length(divenum_CTDtags2009) % ashley is there a way to have 2 figures going at once? or just do 2 loops wie below.
    divenuml = divenum == divenum_CTDtags2009(i);
    ToPlot = matrixCTD(divenuml,:);
    
    oxy = ToPlot(:,7); % oxygen
    dep = ToPlot(:,5); % depth
    la = ToPlot(1,3); % lat
    lo = ToPlot(1,4); % lon
    dat = ToPlot(1,2); % date
    maxdep = max(dep);
    maxdepind = find(maxdep == dep);
    lavector = [lavector; la];
    lovector = [lovector; lo];
    datvector = [datvector; dat];
    
    x = oxy(maxdepind:end);
    y = dep(maxdepind:end);
    oxyvector = [oxyvector; x];
    depvector = [depvector; y];
    
    plot(x,y, 'LineWidth', 4, 'Color', colorCTD(i,:));
    hold on
    text(maxo*0.6, maxd*0.6+(incrd*i), datestr(dat,1), 'Color', colorCTD(i,:),'FontSize', 14, 'FontWeight', 'bold')
end
plot(repmat(.5,1,maxd), 1:maxd, '--', 'LineWidth', 2)
text(maxo*0.6, maxd*0.9, ['n = ' num2str(length(divenum_CTDtags2009)) ' casts '], ...
    'Color', 'k','FontSize', 14, 'FontWeight', 'bold')
hold off
set(gca, 'YDir', 'reverse')
set(gca, 'fontsize', 12, 'fontweight', 'bold')
xlabel('Oxygen (ml/L)','FontSize', 14, 'FontWeight', 'bold')
ylabel('Depth (m)','FontSize', 14, 'FontWeight', 'bold')
title([castID ' CTD Oxygen Comparisons for 2009 Tags, Monterey Bay '], 'FontSize', 18, 'FontWeight', 'bold');
FigName = [castID 'CTDplot_tags2009_Oxy.pdf'];
annotate_JS(Mfilename, gcf, FigName);
%     orient landscape % save
%     print('-dpdf', [dirCTD FigName]);

%% stats

b_ROVCTD_stats % by J. Stewart
%     orient landscape % save
%     print('-dpdf', [dirCTD FigName]);

%% map with locations and dates in coordinating color

figure
mapSqArea
set(gca, 'XLim', lonx, 'YLim', laty);
incrl = range(lonx)*0.02; %this is an el not a one
for i = 1:length(divenum_CTDtags2009)
    hold on
    plot(lovector(i), lavector(i), '.', 'MarkerSize', 20, 'Color', colorCTD(i,:))
    text(lonx(1)+(0.5*range(lonx)), laty(2)-(0.2*range(laty))-(incrl*i), datestr(datvector(i),1), 'Color', colorCTD(i,:), ...
        'FontSize', 14, 'FontWeight', 'bold')
end
plot(LocTagDeployCA(:,1), LocTagDeployCA(:,2), 'r*', 'MarkerSize', 12)
plot(LocTagPopUpCA(:,1), LocTagPopUpCA(:,2), 'r*', 'MarkerSize', 12)
text(lonx(1)+(0.7*range(lonx)), laty(2)-(0.9*range(laty)), '* tag deployments and popoffs' , 'Color', 'r', ...
    'FontSize', 11, 'FontWeight', 'bold')
text(lonx(1)+(0.7*range(lonx)), laty(2)-(0.95*range(laty)), ['n = ' num2str(length(divenum_CTDtags2009)) ' casts '], ...
    'Color', 'k','FontSize', 14, 'FontWeight', 'bold')
set(gca, 'fontsize', 12, 'fontweight', 'bold')
xlabel('Longitude','FontSize', 14, 'FontWeight', 'bold')
ylabel('Latitude','FontSize', 14, 'FontWeight', 'bold')
title([castID ' CTD Locations for 2009 Tags, Monterey Bay'], 'FontSize', 18, 'FontWeight', 'bold');
FigName = [castID 'CTDplot_tags2009_Map.pdf'];
annotate_JS(Mfilename, gcf, FigName);
% orient landscape % save
% print('-dpdf', [dirCTD FigName]);

%% compare temperatures across CTD casts for 2009

figure
maxt = 20; % maxd defined above
for i = 1:length(divenum_CTDtags2009)
    divenuml = divenum == divenum_CTDtags2009(i);
    ToPlot = matrixCTD(divenuml,:);
    
    oxy = ToPlot(:,7); % oxygen
    dep = ToPlot(:,5); % depth
    tem = ToPlot(:,6); % temp
    la = ToPlot(1,3); % lat
    lo = ToPlot(1,4); % lon
    dat = ToPlot(1,2); % date
    maxdep = max(dep);
    maxdepind = find(maxdep == dep);
    
    x = tem(maxdepind:end);
    y = dep(maxdepind:end);
    plot(x,y, 'LineWidth', 4, 'Color', colorCTD(i,:));
    text(maxt*0.6, maxd*0.6+(incrd*i), datestr(dat,1), 'Color', colorCTD(i,:),'FontSize', 14, 'FontWeight', 'bold')
    hold on
end
text(maxt*0.6, maxd*0.9, ['n = ' num2str(length(divenum_CTDtags2009)) ' casts '], ...
    'Color', 'k','FontSize', 14, 'FontWeight', 'bold')
set(gca, 'YDir', 'reverse')%, 'YLim', [0 maxd+50]);
xlabel('Temperature (°C)','FontSize', 14, 'FontWeight', 'bold')
ylabel('Depth (m)','FontSize', 14, 'FontWeight', 'bold')
set(gca, 'fontsize', 12, 'fontweight', 'bold')
if TagTempOverlay % this defined above
    plot(TempTag, DepthTag, '.', 'Color', [0 0.5020 0]); % this line is from processSeriesData.m % J. Stewart
    text(maxt*0.6, maxd*0.9, ['* temp data, tag #' tagNumS '   '], 'Color', [0 0.5020 0], 'FontSize', 11, 'FontWeight', 'bold')
    FigName = [castID 'CTDplot_tags2009_Temp_' tagNumS '.pdf'];
else
    FigName = [castID 'CTDplot_tags2009_Temp.pdf'];
end
title([castID ' CTD Temp Comparisons for 2009 Tags, Monterey Bay  '], 'FontSize', 18, 'FontWeight', 'bold');
annotate_JS(Mfilename, gcf, FigName);
%     orient landscape % save
%     print('-dpdf', [dirCTD FigName]);

% map is same as above from oxygen: don't make a second one.

%% compare temperatures to oxygen across CTD casts for 2009

figure
for i = 1:length(divenum_CTDtags2009)
    divenuml = divenum == divenum_CTDtags2009(i);
    ToPlot = matrixCTD(divenuml,:);
    
    oxy = ToPlot(:,7); % oxygen
    dep = ToPlot(:,5); % depth
    tem = ToPlot(:,6); % temp
    la = ToPlot(1,3); % lat
    lo = ToPlot(1,4); % lon
    dat = ToPlot(1,2); % date
    maxdep = max(dep);
    maxdepind = find(maxdep == dep);
    
    x = oxy(maxdepind:end);
    y = tem(maxdepind:end);
    plot(x,y, 'LineWidth', 4, 'Color', colorCTD(i,:));
    text(maxo*0.6, maxt*0.4-(incrt*i), datestr(dat,1), 'Color', colorCTD(i,:),'FontSize', 14, 'FontWeight', 'bold')
    hold on
end
plot(repmat(.5,1,18), 1:18, '--', 'LineWidth', 2)
text(maxo*0.6, maxt*0.1, ['n = ' num2str(length(divenum_CTDtags2009)) ' casts '], ...
    'Color', 'k','FontSize', 14, 'FontWeight', 'bold')
xlabel('Oxygen (ml/L)','FontSize', 14, 'FontWeight', 'bold')
ylabel('Temperature (°C)','FontSize', 14, 'FontWeight', 'bold')
set(gca, 'fontsize', 12, 'fontweight', 'bold')

title([castID ' CTD Temp-Oxygen Comparisons for 2009 Tags, Monterey Bay  '], 'FontSize', 18, 'FontWeight', 'bold');
FigName = [castID 'CTDplot_tags2009_TempOxy.pdf'];
annotate_JS(Mfilename, gcf, FigName);
% orient landscape % save
% print('-dpdf', [dirCTD FigName]);

%%

disp('Completed plotCTD_Tags_2009.m')
% ===== EOF [plotCTD_Tags_2009.m] ======
