% plotCTD_Tags_Area.m
%
% Called from manageTagData_CTD
%
% Section Titles:
%         Determine which CTDs to use, based on area/location
%         plot oxygen of CTDs in the tag area
%         compare temperatures across CTD casts for area
%
% Outside Functions Called:
% mapSqArea by J. Stewart: maps Monterey area of 2009 tags
% annotate_JS by J. Stewart
% b_ROVCTD_stats % by J. Stewart
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 01-Feb-2011 18:23:22
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : plotCTD_Tags_Area.m

Mfilename = mfilename;

plotspec = 0; % make this user specified, so that plotspec =1 means to specify the month(s) of plotting. So could do upwelling seaon, for ex.
monthspec = 11;
plotfigs = 0;

%% more setup

% plot(SpatialLon, SpatialLat) % it's a box around the tagged area

% identify CTDs in this zone: will have to do a second screening process below because this will include ROV dives that don't start within this area
lat_TagArealog = latCTD > min(SpatialLat) & latCTD < max(SpatialLat);
lon_TagArealog = lonCTD > min(SpatialLon) & lonCTD < max(SpatialLon);

TagArea_Temp = lat_TagArealog + lon_TagArealog;
TagArea_Log = TagArea_Temp == 2;

% redefined matrixCTD that is now just for the specified area
matrixCTD_Area = matrixCTD(TagArea_Log,:); % maybe this isn't even necessary?

divenumtemp = divenum(TagArea_Log);
divenumDA = unique(divenumtemp); % these are the divenums to plot in the area.

runthrough = 0;
dateYrUni = unique(dateCTDstartV(:,1));

%% no plotting here, just go through and make matrices/vectors

lat_vector = [];
lon_vector = [];
dnum_vector = [];
dat_vector = [];
oxyvector = [];
depthvector = [];
temvector = [];
divenumD = [];
for i = 1:length(divenumDA)
    runthrough = runthrough + 1;
    divenuml = divenum == divenumDA(i);
    ToPlot = matrixCTD(divenuml,:);
    
    % need another check since we are indexing ROVs that started inside this zone:
    la = ToPlot(1,3); % lat
    lo = ToPlot(1,4); % lon
    if la > min(SpatialLat) && la < max(SpatialLat) ...
            && lo > min(SpatialLon) && lo < max(SpatialLon);
        oxy = ToPlot(:,7); % oxygen
        dep = ToPlot(:,5); % depth
        tem = ToPlot(:,6);
        dnum = ToPlot(1,1); % divenum
        dat = ToPlot(1,2); % date
        datV = datevec(dat);
        maxdep = max(dep);
        maxdepind = find(maxdep == dep);
        lat_vector = [lat_vector; la];
        lon_vector = [lon_vector; lo];
        dnum_vector = [dnum_vector; dnum];
        dat_vector = [dat_vector; dat]; % this spans from 1999-2009. All months rep'd
        divenumD = [divenumD divenumDA(i)];
        
        x = oxy; % oxy(maxdepind:end);
        y = dep; % dep(maxdepind:end);
        oxyvector = [oxyvector; x];
        depthvector = [depthvector; y];
        temvector = [temvector; tem];
    end
end

dateCTDV = datevec(dateCTD);
dat_vectorV = datevec(dat_vector);% this spans from 1999-2009. All months rep'd

% save -mat rovStats.mat oxyvector depthvector temvector

%%

b_ROVCTD_stats % by J. Stewart; this will also call plotCTD_TagsOverlay

%% figure out spatial component for color
    
    if colorData == 1 % color by years
        colorCTD = lines(length(dateYrUni));
    elseif colorData == 2 % color gradient from shore
        colorincr = 6; % need to divide the lon into colorincr strips and color whatever is inside the strips
        colorstrip = range(lon_vector)/colorincr;
        colorstriptrix = [];
        for i = 1:colorincr
            colorstriptrix = [colorstriptrix min(lon_vector)+i*colorstrip];
        end
        colorstriptrix = [min(lon_vector) colorstriptrix];
        colorCTD = flipud(winter(length(divenumD))); %$$$$$$$$$ the logic is working, but need to make the colors a bit easier to see...
        % also need to change the settings a bit on the plots (although not a
        % huge deal). Make this more robust.
    end
    
%% map with locations and dates in coordinating color
    
    Mfilename = mfilename;
    
    figure
    mapSqArea
    set(gca, 'XLim', lonx, 'YLim', laty);
    incrl = range(lonx)*0.01; %this is an el not a one
    for i = 1:length(divenumD)
        hold on
%         if colorData == 1
%             j = find(dat_vectorV(i,1) == dat_vectorYrs);
%         elseif colorData == 2
%             jtemp = find(lon_vector(i,1) <= colorstriptrix);
%             j = jtemp(1);
%         end
        plot(lon_vector(i), lat_vector(i), '.', 'MarkerSize', 20)%, 'Color', colorCTD(j,:));
        %     text(lon_vector(i), lat_vector(i), num2str(dnum_vector(i)));
    end
    plot(LocTagDeployCA(:,1), LocTagDeployCA(:,2), 'r*', 'MarkerSize', 12)
    plot(LocTagPopUpCA(:,1), LocTagPopUpCA(:,2), 'b*', 'MarkerSize', 12)
    if Spatial == 1
        plot(SpatialLon, SpatialLat) % add the area box around it
    end
%     for j = 1:length(dat_vectorYrs)
%         text(lonx(2)-(0.3*range(lonx)), laty(2)-(0.05*range(laty))-(incrl*j),... %num2str(dat_vectorYrs(j)), ...
%             'Color', colorCTD(j,:), 'FontSize', 14, 'FontWeight', 'bold')
%     end
    text(lonx(2)-(0.3*range(lonx)), laty(2)-(0.9*range(laty)), '* tag deployments and popoffs' , 'Color', 'r', ...
        'FontSize', 11, 'FontWeight', 'bold')
    text(lonx(2)-(0.3*range(lonx)), laty(2)-(0.94*range(laty)), ['n = ' num2str(length(divenumD)) ' casts '], ...
        'Color', 'k','FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'fontsize', 12, 'fontweight', 'bold')
    xlabel('Longitude','FontSize', 14, 'FontWeight', 'bold')
    ylabel('Latitude','FontSize', 14, 'FontWeight', 'bold')
    title([castID ' CTD Locations for Tag Area ' num2str(Spatial)'], 'FontSize', 18, 'FontWeight', 'bold');
    FigName = [castID 'CTDplot_TagArea' num2str(Spatial) '_Map.pdf'];
    annotate_JS(Mfilename, gcf, FigName);
    orient landscape % save
    print('-dpdf', [dirCTD FigName]);
    
    %% To look specifically at upwelling season
if plotfigs    
    uptime = 3:8; % upwelling months
    downtime = [1:2, 9:12]; % davidson (downwelling) months
    dat_vectorV = datevec(dat_vector);
    uptimelog = dat_vectorV(:,2) == uptime(1) | dat_vectorV(:,2) == uptime(2)| dat_vectorV(:,2) == uptime(3) | ...
        dat_vectorV(:,2) == uptime(4) | dat_vectorV(:,2) == uptime(5) | dat_vectorV(:,2) == uptime(6);
    downtimelog = dat_vectorV(:,2) == downtime(1) | dat_vectorV(:,2) == downtime(2)| dat_vectorV(:,2) == downtime(3) | ...
        dat_vectorV(:,2) == downtime(4) | dat_vectorV(:,2) == downtime(5) | dat_vectorV(:,2) == downtime(6);
    
    divenumDUp = divenumD(uptimelog);
    divenumDDown = divenumD(downtimelog);
    
    % Then search: divenum D and replace with divenum DUp or divenum DDown and change the
    % titles. Be sure to also switch back the two lines above involving
    % uptimelog/downtimelog.
    
    %% plot oxygen of CTDs in the tag area
    
    figure
    for i = 1:length(divenumD)
        runthrough = runthrough + 1;
        divenuml = divenum == divenumD(i);
        ToPlot = matrixCTD(divenuml,:);
        
        
        % need another check since we are indexing ROVs that started inside this zone:
        la = ToPlot(1,3); % lat
        lo = ToPlot(1,4); % lon
        if la > min(SpatialLat) && la < max(SpatialLat) ...
                && lo > min(SpatialLon) && lo < max(SpatialLon);
            oxy = ToPlot(:,7); % oxygen
            dep = ToPlot(:,5); % depth
            dat = ToPlot(1,2); % date
            datV = datevec(dat);
            maxdep = max(dep);
            maxdepind = find(maxdep == dep);
            
            x = oxy; % oxy(maxdepind:end);
            y = dep; % dep(maxdepind:end);
            
            colorCTDyrsOrder = [];
            if runthrough > 1
                j = find(datV(1) == dateYrUni); % find index number to identify color
                plot(x,y, 'LineWidth', 4, 'Color', colorCTD(j,:));
                colorCTDyrsOrder = [colorCTDyrsOrder; colorCTD(j,:)];
            else
                plot(x,y, 'LineWidth', 4)
            end
            hold on
        end
    end
    plot(repmat(.5,1,maxd), 1:maxd, '--k', 'LineWidth', 2)
    hold off
    % figure out colors for date vector
    dat_vectorV = datevec(dat_vector);
    dat_vectorYrs = unique(dat_vectorV(:,1));
    for j = 1:length(dat_vectorYrs)
        jdex = find(dat_vectorYrs(j) == dateYrUni);
        text(maxo*0.9, ((maxd+1000)*0.8)+((2*incrd)*j), num2str(dat_vectorYrs(j)), 'Color',...
            colorCTD(jdex,:), 'FontSize', 14, 'FontWeight', 'bold')
    end
    text(maxo*0.9, maxd+1800, ['n = ' num2str(length(divenumD)) ' casts '], ...
        'Color', 'k','FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'YDir', 'reverse')%, 'YLim', [0 maxd+50]);
    set(gca, 'fontsize', 12, 'fontweight', 'bold')
    xlabel('Oxygen (ml/L)','FontSize', 14, 'FontWeight', 'bold')
    ylabel('Depth (m)','FontSize', 14, 'FontWeight', 'bold')
    title([castID ' CTD Oxygen Comparisons for Tag Area ' num2str(Spatial)'], 'FontSize', 18, 'FontWeight', 'bold');
    FigName = [castID 'CTDplot_TagArea_Oxy' num2str(Spatial) '.pdf'];
    annotate_JS(Mfilename, gcf, FigName);
    % orient landscape % save
    % print('-dpdf', [dirCTD FigName]);
    
    % plot oxygen of CTDs in the tag area 
    
    
    % weird number on map figure
    % figure
    % for i = 1:length(divenumD)
    %     hold on
    %     plot(lon_vector(i), lat_vector(i), '.', 'MarkerSize', 1, 'Color', colorCTD(j,:));
    %     text(lon_vector(i), lat_vector(i), num2str(dnum_vector(i)));
    % end
    
    %% compare temperatures across CTD casts for area
    
    figure
    plotCount = 0;
    maxt = 20; % maxd defined above
    for i = 1:length(divenumD)
        divenuml = divenum == divenumD(i);
        ToPlot = matrixCTD(divenuml,:);
        
        dat = ToPlot(1,2); % date
        datV = datevec(dat);
        if plotspec
            if datV(:,2) == monthspec
                plotCount = plotCount + 1;
                
                oxy = ToPlot(:,7); % oxygen
                dep = ToPlot(:,5); % depth
                tem = ToPlot(:,6); % temp
                la = ToPlot(1,3); % lat
                lo = ToPlot(1,4); % lon
                
                datV = datevec(dat);
                maxdep = max(dep);
                maxdepind = find(maxdep == dep);
                
                x = tem;
                y = dep;
                j = find(datV(1) == dat_vectorYrs);
                plot(x,y, 'LineWidth', 4, 'Color', colorCTD(j,:));
                hold on
            end
        else
            plotCount = length(divenumD);
            oxy = ToPlot(:,7); % oxygen
            dep = ToPlot(:,5); % depth
            tem = ToPlot(:,6); % temp
            la = ToPlot(1,3); % lat
            lo = ToPlot(1,4); % lon
            
            datV = datevec(dat);
            maxdep = max(dep);
            maxdepind = find(maxdep == dep);
            
            x = tem;
            y = dep;
            j = find(datV(1) == dat_vectorYrs);
            plot(x,y, 'LineWidth', 4, 'Color', colorCTD(j,:));
            hold on
        end
    end
    for j = 1:length(dat_vectorYrs)
        text(maxt*0.8, ((maxd+1000)*0.8)+((2*incrd)*j), num2str(dat_vectorYrs(j)), ...
            'Color', colorCTD(j,:), 'FontSize', 14, 'FontWeight', 'bold')
    end
    grid on
    text(maxt*0.8, maxd+1800, ['n = ' num2str(plotCount) ' casts '], ...
        'Color', 'k','FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'YDir', 'reverse')%, 'YLim', [0 maxd+50]);
    xlabel('Temperature (°C)','FontSize', 14, 'FontWeight', 'bold')
    ylabel('Depth (m)','FontSize', 14, 'FontWeight', 'bold')
    set(gca, 'fontsize', 12, 'fontweight', 'bold')
    if TagTempOverlay % this defined above
        plot(TempTag, DepthTag, '.', 'Color', 'k'); % this line is from processSeriesData.m % J. Stewart
        text(maxt*0.6, (maxd*0.9)+2000, ['* temp data, tag #' tagNumS '   '], 'Color', 'k', 'FontSize', 11, 'FontWeight', 'bold')
        FigName = [castID 'CTDplot_TagArea_Temp_' tagNumS '' num2str(Spatial) '.pdf'];
    else
        FigName = [castID 'CTDplot_TagArea_Temp' num2str(Spatial) '.pdf'];
    end
    title([castID ' CTD Temp Comparisons for Tag Area ' num2str(Spatial)], 'FontSize', 18, 'FontWeight', 'bold');
    annotate_JS(Mfilename, gcf, FigName);
    % orient landscape % save
    % print('-dpdf', [dirCTD FigName]);
    
    % troubleshoot
    % klog = matrixCTD(:,5) < 50 & matrixCTD(:,6) < 5;
    % sum(klog)
    % a = matrixCTD(klog,1);
    % datestr(matrixCTD(klog,2))
    % figure; hold on
    % for i = 1:length(unique(a))
    %     ilog = a(i) == matrixCTD(:,1);
    %     plot(matrixCTD(ilog,6),matrixCTD(ilog,5), '.')
    % end
    
    %% compare temperatures to oxygen across CTD casts across area
    
    figure
    for i = 1:length(divenumD)
        divenuml = divenum == divenumD(i);
        ToPlot = matrixCTD(divenuml,:);
        
        oxy = ToPlot(:,7); % oxygen
        dep = ToPlot(:,5); % depth
        tem = ToPlot(:,6); % temp
        la = ToPlot(1,3); % lat
        lo = ToPlot(1,4); % lon
        dat = ToPlot(1,2); % date
        datV = datevec(dat);
        maxdep = max(dep);
        maxdepind = find(maxdep == dep);
        
        x = oxy;
        y = tem;
        j = find(datV(1) == dat_vectorYrs);
        plot(x,y, '.', 'MarkerSize', 15, 'Color', colorCTD(j,:));
        hold on
    end
    for j = 1:length(dat_vectorYrs)
        text(maxo*0.9, maxt*0.4-(incrt*j), num2str(dat_vectorYrs(j)), ...
            'Color', colorCTD(j,:), 'FontSize', 14, 'FontWeight', 'bold')
    end
    text(maxo*0.9, maxt*0.05, ['n = ' num2str(length(divenumD)) ' casts '], ...
        'Color', 'k','FontSize', 14, 'FontWeight', 'bold');
    plot(repmat(.5,1,18), 1:18, '--', 'LineWidth', 2)
    xlabel('Oxygen (ml/L)','FontSize', 14, 'FontWeight', 'bold')
    ylabel('Temperature (°C)','FontSize', 14, 'FontWeight', 'bold')
    set(gca, 'fontsize', 12, 'fontweight', 'bold')
    title([castID ' CTD Temp-Oxygen Comparisons for Tag Area ' num2str(Spatial)], 'FontSize', 18, 'FontWeight', 'bold');
    FigName = [castID 'CTDplot_TagArea_TempOxy' num2str(Spatial) '.pdf'];
    annotate_JS(Mfilename, gcf, FigName);
    % orient landscape % save
    % print('-dpdf', [dirCTD FigName]);
    
end
%%

disp('Completed plotCTD_Tags_Area.m')
% ===== EOF [plotCTD_Tags_Area.m] ======
