% importCalCOFIbottles.m
%
% Called from manageTagData_CTD.m
%
% Outside Functions Called:
%     sw_dpth.m from CISRO http://www.cmar.csiro.au/datacentre/ext_docs/seawater.htm
%       plotCalCOFI_lines.m by J. Stewart
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 06-Apr-2010 16:50:41
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : importCalCOFIbottles.m

% to read in CalCOFI water bottle samples
% CalCOFI cruise went from south to north (grah!)
%identify station unique to each date
% line 83 station 55,60 both on 11/17/09 (=day13) so that's a few days after pass Pt Conception

% first station 49 is N of pt conception

Mfilename = mfilename;
userinput = 0;
figs = 1;
TagTempOverlay = 0;

% incorporate this from email (search bograd isotherm): two text files with columns as following: lon, lat,
% depth of 10C isotherm, line, station for the long-term mean (84-09) and Nov 2009 cruise individually.

%% Read in File

cd '~/Dropbox/Thesis/Dg_Tagging/_CTDs/_CalCOFI_CTD_data/';
filename = '~/Dropbox/Thesis/Dg_Tagging/_CTDs/_CalCOFI_CTD_data/calcofi0911_final.txt';
castID = 'CalCOFIbottle';

fid=fopen(filename);
line = 1;
% while line<5 %get headers
tline = fgetl(fid);
if ~isempty(tline) && strcmp(tline(1),'D'), break, end %get header names
%  line = line+1;
%  end
fclose(fid);
Tline = ['''' regexprep(tline,'\t',''';''') '''']; %reformat headers
CalCOFI.headernames = eval(['{' Tline '}']); %make into string array
CalCOFI.headernames = CalCOFI.headernames';
numcol = length(CalCOFI.headernames);

%import column data
fmat = '%n'; %initial format of each columm
for s = 1:numcol %make colheaderFormat line adding on %n as long as headers are
    fmat = strcat(fmat,' %n');
end
fid=fopen(filename);
CalCOFI.data = textscan(fid, fmat,'HeaderLines', line, 'Delimiter','\t');
fclose(fid);

clear tline fid numcol numcols Tline

%Create CastNumberID, wie divenum
% for i = 1:size(UniCast,1)
%     linlog = UniCast(i,1) == LineNum;
%     stalog = UniCast(i,2) == StationNum;
%     linstalog = linlog+stalog;
%     linsta = linstalog == 2;
%
% end
% CastNumberID % equivalent of divenum

%% Setup

% redefine
DateCal = datenum(CalCOFI.data{1,1:5},0); % convert dates
LatCal = CalCOFI.data{1,8};
LonCal = -CalCOFI.data{1,9};
LineNum = CalCOFI.data{1,6};
StationNum = CalCOFI.data{1,7};
OxyCal_mlL = CalCOFI.data{1,16}; OxyCal_mlL(OxyCal_mlL == -999) = NaN;
TempCal_C = CalCOFI.data{1,12};
DepthCal_m = sw_dpth(CalCOFI.data{1,11}, LatCal); % this had been pressure

% ID casts
UniLineNum = unique(LineNum); % unique CalCOFI line numbers
UniStationNum = unique(StationNum); % unique CalCOFI station numbers
CastAll = [LineNum StationNum];
UniCast = unique(CastAll, 'rows');
Coords = [LonCal LatCal];
UniCoords = unique(Coords, 'rows');

% DataCleanUp % strange data >2000m, but this will get rid of entire lines and cause problems.
% a = find(DepthCal_m > 2000);
% problemCasts = CastAll(a,:); %80-90 and 90-90
% problemCastsLine = unique(problemCasts(:,1));
% problemCastsStation = unique(problemCasts(:,2));
% problemCastsLogA = LineNum == problemCastsLine(1) | LineNum == problemCastsLine(2); % if any more casts, will have to account for that
% problemCastsLogB = StationNum == problemCastsStation;
% problemCastsLogAB = problemCastsLogA + problemCastsLogA;
% problemCastsLog = problemCastsLogAB == 2;
% DateCal(problemCastsLog) = []; % cleanup
% LatCal(problemCastsLog) = [];
% LonCal(problemCastsLog) = [];
% LineNum(problemCastsLog) = [];
% StationNum(problemCastsLog) = [];
% OxyCal_mlL(problemCastsLog) = [];
% TempCal_C(problemCastsLog) = [];
% DepthCal_m(problemCastsLog) = [];
% UniLineNum = unique(LineNum); % unique CalCOFI line numbers % yes, again.
% UniStationNum = unique(StationNum); % unique CalCOFI station numbers
% CastAll = [LineNum StationNum];
% UniCast = unique(CastAll, 'rows');

% colormap
colorCalCOFIa = flipud(winter(length(UniStationNum)));
colorCalCOFI = colorCalCOFIa(1:length(UniStationNum),:);

% save to import into manageCTD_Tags
% matrixCalCOFICTD = [CalCOFICTD.data{1,1} DateCal LatCal LonCal DepthCal_m TempCal_C...
%       OxyCal_mlL];
% matrixCalCOFICTDheaders = {'divenum' 'dateCTD' 'lat' 'lon' 'depth' 'temper' 'oxy_mlL'};
%
% cd('/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/_CTDs/_CalCOFI_CTD_data');
% save -mat CalCOFICTD.mat matrixCalCOFICTD matrixCalCOFICTDheaders

% save all this in a structure
for i = 1:length(UniLineNum)
    ilog = LineNum == UniLineNum(i);
    for j = 1:length(UniStationNum) % for each station in that line
        jlog = StationNum == UniStationNum(j);
        ij = ilog + jlog;
        ijlog = ij == 2;
        
        CalCOFIcasts(i,j).date = DateCal(ijlog);
        CalCOFIcasts(i,j).temp = TempCal_C(ijlog);
        CalCOFIcasts(i,j).depth = DepthCal_m(ijlog);
        CalCOFIcasts(i,j).line = UniLineNum(i);
        CalCOFIcasts(i,j).station = UniStationNum(j);
        CalCOFIcasts(i,j).lon = LonCal(ijlog);
        CalCOFIcasts(i,j).lat = LatCal(ijlog);
    end
end


%% plot all temp, all oxy, temp v. oxy, then map.

if figs
    
    maxt = 20;
    maxd = 800;
    maxo = 8;
    incrt = maxt * 0.5; % temp increment
    incrd = maxd * 0.05; % temp increment
    
    % temper
    figure
    hold on
    for i = 1:length(UniStationNum)
        ilog = UniStationNum(i) == StationNum;
        plot(TempCal_C(ilog), DepthCal_m(ilog), '.', 'Color', colorCalCOFI(i,:), 'MarkerSize', 15);
    end
    set(gca, 'YDir', 'reverse', 'YLim', [0 maxd]);
    text(maxt*0.4, maxd*0.8, ['n = ' num2str(size(UniCast,1)) ' casts '], ...
        'Color', 'k','FontSize', 14, 'FontWeight', 'bold');
    text(maxt*0.4, maxd*0.85, 'color gradient: land (green) to offshore (blue)', ...
        'Color', 'k','FontSize', 14, 'FontWeight', 'bold');
    xlabel('Temperature (°C)','FontSize', 14, 'FontWeight', 'bold')
    ylabel('Depth (m)','FontSize', 14, 'FontWeight', 'bold')
    set(gca, 'fontsize', 12, 'fontweight', 'bold')
    if TagTempOverlay % this defined above
        plot(TempTag, DepthTag, '.', 'Color', 'k'); % this line is from processSeriesData.m % J. Stewart
        text(maxt*0.6, (maxd*0.9)+2000, ['* temp data, tag #' tagNumS '   '], 'Color', 'k', 'FontSize', 11, 'FontWeight', 'bold')
        FigName = [castID 'CTDplot_TagArea_Temp_' tagNumS '.pdf'];
    else
        FigName = [castID 'CTDplot_TagArea_Temp.pdf'];
    end
    title([castID ' CTD Temp Comparisons '], 'FontSize', 18, 'FontWeight', 'bold');
    annotate_JS(Mfilename, gcf, FigName);
    % orient landscape % save
    % print('-dpdf', [dirCTD FigName]);
    
    % oxy
    figure
    hold on
    for i = 1:length(UniStationNum)
        ilog = UniStationNum(i) == StationNum;
        plot(OxyCal_mlL(ilog), DepthCal_m(ilog), '.', 'Color', colorCalCOFI(i,:), 'MarkerSize', 15);
    end
    set(gca, 'YDir', 'reverse', 'YLim', [0 maxd]);
    text(maxo*0.2, maxd*0.8, ['n = ' num2str(size(UniCast,1)) ' casts '], ...
        'Color', 'k','FontSize', 14, 'FontWeight', 'bold');
    text(maxo*0.2, maxd*0.85, 'color gradient: land (green) to offshore (blue)', ...
        'Color', 'k','FontSize', 14, 'FontWeight', 'bold');
    plot(repmat(.5,1,maxd), 1:maxd, '--k', 'LineWidth', 2)
    set(gca, 'YDir', 'reverse', 'YLim', [0 maxd]);
    xlabel('Oxygen (mL/L)','FontSize', 14, 'FontWeight', 'bold')
    ylabel('Depth (m)','FontSize', 14, 'FontWeight', 'bold')
    set(gca, 'fontsize', 12, 'fontweight', 'bold')
    if TagTempOverlay % this defined above
        plot(TempTag, DepthTag, '.', 'Color', 'k'); % this line is from processSeriesData.m % J. Stewart
        text(maxt*0.6, (maxd*0.9)+2000, ['* temp data, tag #' tagNumS '   '], 'Color', 'k', 'FontSize', 11, 'FontWeight', 'bold')
        FigName = [castID 'CTDplot_TagArea_Oxy_' tagNumS '.pdf'];
    else
        FigName = [castID 'CTDplot_TagArea_Oxy.pdf'];
    end
    title([castID ' CTD Oxygen Comparisons '], 'FontSize', 18, 'FontWeight', 'bold');
    annotate_JS(Mfilename, gcf, FigName);
    % orient landscape % save
    % print('-dpdf', [dirCTD FigName]);
    
    % temp-oxy
    figure
    hold on
    for i = 1:length(UniStationNum)
        ilog = UniStationNum(i) == StationNum;
        plot(OxyCal_mlL(ilog), TempCal_C(ilog), '.', 'Color', colorCalCOFI(i,:), 'MarkerSize', 15);
    end
    set(gca, 'YLim', [0 maxt]);
    text(maxo*0.2, maxt*0.15, ['n = ' num2str(size(UniCast,1)) ' casts '], ...
        'Color', 'k','FontSize', 14, 'FontWeight', 'bold');
    text(maxo*0.2, maxt*0.2, 'color gradient: land (green) to offshore (blue)', ...
        'Color', 'k','FontSize', 14, 'FontWeight', 'bold');
    plot(repmat(.5,1,maxt), 1:maxt, '--k', 'LineWidth', 2)
    xlabel('Oxygen (mL/L)','FontSize', 14, 'FontWeight', 'bold')
    ylabel('Temperature (°C)','FontSize', 14, 'FontWeight', 'bold')
    set(gca, 'fontsize', 12, 'fontweight', 'bold')
    if TagTempOverlay % this defined above
        plot(TempTag, DepthTag, '.', 'Color', 'k'); % this line is from processSeriesData.m % J. Stewart
        text(maxt*0.6, (maxd*0.9)+2000, ['* temp data, tag #' tagNumS '   '], 'Color', 'k', 'FontSize', 11, 'FontWeight', 'bold')
        FigName = [castID 'CTDplot_TagArea_TempOxy_' tagNumS '.pdf'];
    else
        FigName = [castID 'CTDplot_TagArea_TempOxy.pdf'];
    end
    title([castID ' CTD Temp-Oxy Comparisons '], 'FontSize', 18, 'FontWeight', 'bold');
    annotate_JS(Mfilename, gcf, FigName);
    % orient landscape % save
    % print('-dpdf', [dirCTD FigName]);
    
    %% map it please
    
    % SpatialLon = [-117; -124.5; -124.5; -117; -117]; % based on min, max LonCal
    % SpatialLat = [29.5; 29.5; 35.5; 35.5; 29.5];
    % Locations = [SpatialLat SpatialLon]; % just defined in plotCTD_Tags_Area
    figure
    mapSqArea
    hold on
    set(gca, 'XLim', [-125 max(SpatialLon)], 'YLim', [min(SpatialLat) 37]);
    % incrl = range(lonx)*0.03; %this is an el not a one
    for i = 1:length(UniStationNum)
        ilog = UniStationNum(i) == StationNum;
        plot(LonCal(ilog), LatCal(ilog), '.', 'Color', colorCalCOFI(i,:), 'MarkerSize', 15);
    end
    plot(LocTagDeployCA(:,1), LocTagDeployCA(:,2), 'r*', 'MarkerSize', 12)
    plot(LocTagPopUpCA(:,1), LocTagPopUpCA(:,2), 'k*', 'MarkerSize', 12)
    text(min(SpatialLon)+(0.3*range(SpatialLon)), min(SpatialLat)+(0.1*range(SpatialLat)), '* tag deployments and popoffs' , 'Color', 'r', ...
        'FontSize', 11, 'FontWeight', 'bold')
    text(min(SpatialLon)+(0.3*range(SpatialLon)), min(SpatialLat)+(0.05*range(SpatialLat)), ['n = ' num2str(size(UniCast,1)) ' casts '], ...
        'Color', 'k','FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'fontsize', 12, 'fontweight', 'bold')
    xlabel('Longitude','FontSize', 14, 'FontWeight', 'bold')
    ylabel('Latitude','FontSize', 14, 'FontWeight', 'bold')
    title([castID ' CTD Locations for CalCOFI and tags'], 'FontSize', 18, 'FontWeight', 'bold');
    FigName = [castID 'CTDplot_TagArea_Map.pdf'];
    annotate_JS(Mfilename, gcf, FigName);
    % orient landscape % save
    % print('-dpdf', [dirCTD FigName]);
    
    
    %% test for missing station num
    figure
    mapSqArea
    hold on
    set(gca, 'XLim', [-125 max(SpatialLon)], 'YLim', [min(SpatialLat)-1 37]);
    plot(UniCoords(1:end-1,1),UniCoords(1:end-1,2), '.m') % 1:end-1 so it's the same length as UniCasts
    for i = 1:length(UniStationNum)
        ilog = UniStationNum(i) == StationNum;
        plot(LonCal(ilog), LatCal(ilog), '.', 'MarkerSize', 15);
        text(LonCal(ilog), LatCal(ilog), num2str(UniStationNum(i)), 'Color', 'k', ...
            'FontSize', 15, 'FontWeight', 'bold')
    end
    
end
%%
disp('Completed importCalCOFIbottles.m')
% ===== EOF [importCalCOFIbottles.m] ======

