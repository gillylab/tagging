% calcofi_tags_ols.m
%
% Make Least Squares Matrices for each each day of the tag (ascents and
% descents)s comparisons to each CalCOFI station
%
% Read in SquidProfileStructure from the appropriate tag, and compare this
% to the CalCOFI grid. End up with values in the structure for the ordinary
% least squares (ols) for each point in the CalCOFI grid (nans where there
% is no CalCOFI profile). Will run and save for ascents and descents.
%
% plotting is down at the bottom
%
% Called From: manageTagData_CTD.m

% Description:
%
% Outside Functions Called:
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 26-Apr-2011 15:31:08
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : calcofi_tags_ols.m

Mfilename = mfilename;
numcols = 3;
numberofpanels = 6;
plothere = 0;
plotOLSprofiles = 1;
b = 2; % buffer

%% import tag and CalCOFI bottle data if you want this stand-alone

if b == 1
    load -mat SquidProfileStructure_83051_1hrbuff.mat
    zz = [18 18; 17 17]; % days that don't work for ascents, descents because of squid behavior
elseif b == 2
    load -mat SquidProfileStructure_83051_2hrbuff.mat % tag that went into CalCOFI grid. created in plotSeriesData.m % J. Stewart
    zz = [9 10; 2 16]; % days that don't work for ascents, descents because of squid behavior
else
    disp('need to load something');
end
% load -mat CalCOFIcasts.mat % CalCOFI bottles

% setup

% this goes through the process of day/night/day/night (descent/ascent/descent/ascent)
LeastSquaresTrixA = nan*zeros(length(UniLineNum),length(UniStationNum));
LeastSquaresTrixD = nan*zeros(length(UniLineNum),length(UniStationNum));
LontrixA = nan*zeros(length(UniLineNum),length(UniStationNum));
LattrixA = nan*zeros(length(UniLineNum),length(UniStationNum));
LontrixD = nan*zeros(length(UniLineNum),length(UniStationNum));
LattrixD = nan*zeros(length(UniLineNum),length(UniStationNum));

nn = 2;

%% ASCENTS

lsvalA = zeros(nn,length(SquidProfile)); % least squares values for each day
lsindA = zeros(nn,length(SquidProfile)); % least squares indices for each day
lsdateA = zeros(nn,length(SquidProfile));

newAtrix = []; % newAtrix = nan(length(UniLineNum)*length(UniStationNum),length(SquidProfile));
% let's do subplot figures with descents and ascents as the 2 suplots. then compare the values.
Count = 0;
panelnum = 0;
for z = 1:length(SquidProfile)
    Count = Count + 1;
    if z ~= zz(1,1) && z ~= zz(1,2) % Ascents: skip profile 9 and 10 - no useable% test j=24, i=1
        if plotOLSprofiles
            fignum = floor((Count - 1)./numberofpanels)+1;
            fignum2 = fignum+60;
            panelnum = rem(Count - 1, numberofpanels)+1; % rem: the remainder after division
            if panelnum == 1
                figure(fignum2); clf;
            end
        end
        for i = 1:length(UniLineNum)
            for j = 1:length(UniStationNum) % for each station in that line
                if ~isempty(CalCOFIcasts(i,j).temp)
                    
                    if max(CalCOFIcasts(i,j).depth) > max(SquidProfile(z).depthAscent) % Ascent
                        
                        % get squid profile and prepare for use
                        Xraw = SquidProfile(z).depthAscent;
                        Yraw = SquidProfile(z).tempAscent;
                        
                        % perform a 1-metre bin average to get a monotonic vector in X (depth)
                        xtmp = floor(min(Xraw)):1:ceil(max(Xraw));
                        ytmp = nan*ones(size(xtmp));
                        for k = 1:length(xtmp),
                            ibin = find(Xraw<=xtmp(k)+.5 & Xraw>xtmp(k)-.5);
                            if ~isempty(ibin)
                                ytmp(k) = nanmean(Yraw(ibin));
                            end
                        end
                        
                        % get rid of empty bins
                        igood = find(~isnan(ytmp));
                        X = xtmp(igood);
                        Y = ytmp(igood);
                        
                        % grab calcofi depths for interpolation
                        Xi = CalCOFIcasts(i,j).depth;
                        
                        % use spline in interp1, and enforce Nan for extrapolated points - note that there should not be any nans at this point, but better safe than sorry.
                        Yi = interp1(X(~isnan(Y) & ~isnan(X)),Y(~isnan(Y) & ~isnan(X)),Xi, 'spline', nan);
                        YY = CalCOFIcasts(i,j).temp;
                        
                        % get rms difference (root-mean-square, also called quadratic mean)
                        dy = Yi(~isnan(Yi)) - YY(~isnan(Yi));
                        a = sqrt(sum(dy.*dy)/length(find(~isnan(Yi)))); %normalized by number of points in each profile
                        LeastSquaresTrixA(i,j) = a;
                        
                        Loni = CalCOFIcasts(i,j).lon;
                        Lati = CalCOFIcasts(i,j).lat;
                        Dati = CalCOFIcasts(i,j).date;
                        LontrixA(i,j) = Loni(1);
                        LattrixA(i,j) = Lati(1);
                        DatetrixA(i,j) = Dati(1);
                        
                        if plotOLSprofiles
                            xx = CalCOFIcasts(i,j).temp;
                            yy = CalCOFIcasts(i,j).depth;
                            
                            subplot(ceil(numberofpanels/numcols), numcols, panelnum); hold on;
                            plot(xx, yy, 'color', [0.9 0.9 0.9]);
                        end
                    end
                end
            end
        end
        SquidProfile(z).LeastSquaresAscent = LeastSquaresTrixA;
        % now pull out indices where they are the lowest (least squares)
        [newA, indA] = sort(LeastSquaresTrixA(:),1, 'ascend');
        %         newAtrix(:,z) = newA;
        newAtrix = [newAtrix; newA];
        lsvalA(:,z) = newA(1:nn);
        lsindA(:,z) = indA(1:nn);
        lsdateA(:,z) = DatetrixA(lsindA(:,z));
        
        if plotOLSprofiles % plot ascents and descents on the same figure
            for g = 1:nn% plot nn closest casts in color % amount of least squares matches you choose
                plot(CalCOFIcasts(lsindA(g,z)).temp, CalCOFIcasts(lsindA(g,z)).depth, 'color', [0.4 0.4 0.4], 'LineWidth', 1.5); % 'color', ccD(g,:)
                text(12, 200+g*40,...
                    [datestr(lsdateA(g,z), 'mmm-dd') ', rms=' num2str(lsvalA(g,z), '%2.2f') ' '], 'color', ccA(g,:))
            end
            title({['Day ' num2str(z-1) ': ' datestr(HinkeDates(z),'mmm-dd') '  ']; }, 'FontSize',10, 'FontWeight', 'bold')
            plot(Y, X, 'k', 'LineWidth', 4)%, 'markersize', 15) % plot smoothed SquidProfile
            set(gca, 'YLim', [0 600], 'YDir', 'reverse')
            if panelnum == numberofpanels
                key = ['CalCOFI (gray) v. Squid ' num2str(TagHinke) ' Ascents (orange)'];
                [h3] = suptitle(key);
                [ax,h(1)] = suplabel(['temp'],'x');
                [ax,h(2)] = suplabel('depth','y');
                set(h3,'FontSize',24, 'FontWeight', 'bold')
                set(h,'FontSize',15, 'FontWeight', 'bold')
                FigName = [num2str(TagHinke) 'olsCalCOFI_A_' num2str(b) 'hrbuff' num2str(fignum) '.pdf'];
                annotate_JS(Mfilename, gcf, FigName)
%                 orient landscape
%                 print('-dpdf', [dirFortran FigName]);
            end
        end
    else
        LeastSquaresTrixA = nan*zeros(length(UniLineNum),length(UniStationNum));
        SquidProfile(z).LeastSquaresAscent = nan*ones(size(LeastSquaresTrixA));
        % now pull out indices where they are the lowest (least squares)
        [newA, indA] = sort(LeastSquaresTrixA(:),1, 'ascend');
%         newAtrix(:,z) = newA;
newAtrix = [newAtrix; newA];
        lsvalA(:,z) = newA(1:nn);
        lsindA(:,z) = indA(1:nn);
        lsdateA(:,z) = DatetrixA(lsindA(:,z));
    end
end
if plotOLSprofiles
    if panelnum < numberofpanels
        key = ['CalCOFI (gray) v. Squid ' num2str(TagHinke) 'Ascents (orange)'];
        [h3] = suptitle(key);
        [ax,h(1)] = suplabel(['temp'],'x');
        [ax,h(2)] = suplabel('depth','y');
        set(h3,'FontSize',24, 'FontWeight', 'bold')
        set(h,'FontSize',15, 'FontWeight', 'bold')
        FigName = [num2str(TagHinke) 'olsCalCOFI_A_' num2str(b) 'hrbuff' num2str(fignum) '.pdf'];
        annotate_JS(Mfilename, gcf, FigName)
%         orient landscape % save
%         print('-dpdf', [dirFortran FigName]);
    end
end

% make lat and lon into a vector with indices instead of matrices
newLonA = [];
newLatA = [];
for k = 1:size(LontrixA,2)
    newLonA = [newLonA; LontrixA(:,k)];
    newLatA = [newLatA; LattrixA(:,k)];
end

% find the date of these casts
newDateA = [];
for k = 1:size(DatetrixA,2)
    newDateA = [newDateA; DatetrixA(:,k)];
end


% plot RMS histos
figure; hold on
[n, xout] = hist(newAtrix,50);
[n2, xout2] = hist(lsvalA(:,11),xout); % for stacking
N = [n2; n; nan(1,length(n2))]'; % to make the colormap work with stacking
bar(N, 'stacked')
colormap(gray)
xlabel('RMS value ','FontSize', 14, 'FontWeight', 'bold')
ylabel('Count', 'FontSize', 14, 'FontWeight', 'bold')
title({[num2str(TagHinke) ' RMS histogram: Best' num2str(nn) 'Ascents ']; 'RMS of All Squid profiles v. all CalCOFI casts'}, 'fontsize', 16, 'fontweight', 'bold');
FigName = [num2str(TagHinke) '_histogram_RMS_A_' num2str(b) 'hrbuff.pdf']; %FigName = [num2str(TagHinke) '_timeseries_distoffshore_perc.pdf'];
annotate_JS(Mfilename, gcf, FigName);
% orient landscape % save
% print('-dpdf', [dirFortran FigName]);

save -mat calcofi_tags_olsA.mat UniCoords SpatialLon SpatialLat lsindA lsvalA newLonA newLatA

%% DESCENTS

lsvalD = zeros(nn,length(SquidProfile)); % least squares values for each day
lsindD = zeros(nn,length(SquidProfile)); % least squares indices for each day
lsdateD = zeros(nn,length(SquidProfile));
newDtrix = [];
Count = 0;
for z = 1:length(SquidProfile)
    Count = Count + 1;
    if z ~= zz(2,1) && z ~= zz(2,2) % Descents: skip profile 2 and 17 - no useable data % test (10,1)
        if plotOLSprofiles
            fignum = floor((Count - 1)./numberofpanels)+1;
            fignum2 = fignum+50;
            panelnum = rem(Count - 1, numberofpanels)+1; % rem: the remainder after division
            if panelnum == 1
                figure(fignum2); clf;
            end
        end
        for i = 1:length(UniLineNum)
            for j = 1:length(UniStationNum) % for each station in that line
                if ~isempty(CalCOFIcasts(i,j).temp)
                    if max(CalCOFIcasts(i,j).depth) > max(SquidProfile(z).depthDescent) % CalCOFI site must be deeper than squid profile
                        
                        % get squid profile and prepare for use
                        Xraw = SquidProfile(z).depthDescent;
                        Yraw = SquidProfile(z).tempDescent;
                        
                        % perform a 1-metre bin average to get a monotonic vector in X (depth)
                        xtmp = floor(min(Xraw)):1:ceil(max(Xraw));
                        ytmp = nan*ones(size(xtmp));
                        for k = 1:length(xtmp),
                            ibin = find(Xraw<=xtmp(k)+.5 & Xraw>xtmp(k)-.5);
                            if ~isempty(ibin)
                                ytmp(k) = nanmean(Yraw(ibin));
                            end
                        end
                        
                        % get rid of empty bins
                        igood = find(~isnan(ytmp));
                        X = xtmp(igood);
                        Y = ytmp(igood);
                        %clear xtmp ytmp
                        
                        % grab calcofi depths for interpolation
                        Xi = CalCOFIcasts(i,j).depth;
                        
                        % use spline in interp1, and enforce Nan for extrapolated points - note that there should not be any nans at this point, but better safe than sorry.
                        Yi = interp1(X(~isnan(Y) & ~isnan(X)),Y(~isnan(Y) & ~isnan(X)),Xi, 'spline', nan);
                        YY = CalCOFIcasts(i,j).temp;
                        
                        % get rms difference (root-mean-square, also called quadratic mean)
                        dy = Yi(~isnan(Yi)) - YY(~isnan(Yi));
                        a = sqrt(sum(dy.*dy)/length(find(~isnan(Yi)))); %normalized by number of points in each profile
                        
                        LeastSquaresTrixD(i,j) = a;
                        
                        Loni = CalCOFIcasts(i,j).lon;
                        Lati = CalCOFIcasts(i,j).lat;
                        Dati = CalCOFIcasts(i,j).date;
                        LontrixD(i,j) = Loni(1);
                        LattrixD(i,j) = Lati(1);
                        DatetrixD(i,j) = Dati(1);
                        
                        if plotOLSprofiles
                            xx = CalCOFIcasts(i,j).temp;
                            yy = CalCOFIcasts(i,j).depth;
                            
                            subplot(ceil(numberofpanels/numcols), numcols, panelnum); hold on;
                            plot(xx, yy, 'color', [0.8 0.8 0.8]);
                        end
                    end
                end
            end
        end
        SquidProfile(z).LeastSquaresDescent = LeastSquaresTrixD;
        % now pull out indices where they are the lowest (least squares)
        [newD, indD] = sort(LeastSquaresTrixD(:), 1, 'ascend');
        newDtrix = [newDtrix; newD];
        lsvalD(:,z) = newD(1:nn);
        lsindD(:,z) = indD(1:nn);
        lsdateD(:,z) = DatetrixD(lsindD(:,z));
        
        if plotOLSprofiles % plot ascents and descents on the same figure
            for g = 1:nn% plot nn closest casts in color % amount of least squares matches you choose
                plot(CalCOFIcasts(lsindD(g,z)).temp, CalCOFIcasts(lsindD(g,z)).depth, 'color', [0.4 0.4 0.4], 'LineWidth', 1.5); % 'color', ccD(g,:)
                text(12, 400+g*40,...
                    [datestr(lsdateD(g,z), 'mmm-dd') ', rms=' num2str(lsvalD(g,z), '%2.2f') ' '], 'color', ccD(g,:))
            end
            title({['Day ' num2str(z-1) ': ' datestr(HinkeDates(z),'mmm-dd') '  ']; }, 'FontSize',10, 'FontWeight', 'bold')
            plot(Y, X, 'k', 'LineWidth', 4)%, 'markersize', 15) % plot smoothed SquidProfile
            set(gca, 'YLim', [0 600], 'YDir', 'reverse')
            if panelnum == numberofpanels
                key = ['CalCOFI (gray) v. Squid ' num2str(TagHinke) 'Descents (pink)'];
                [h3] = suptitle(key);
                [ax,h(1)] = suplabel(['temp'],'x');
                [ax,h(2)] = suplabel('depth','y');
                set(h3,'FontSize',24, 'FontWeight', 'bold')
                set(h,'FontSize',15, 'FontWeight', 'bold')
                FigName = [num2str(TagHinke) 'olsCalCOFI_D_' num2str(b) 'hrbuff' num2str(fignum) '.pdf'];
                annotate_JS(Mfilename, gcf, FigName)
%                 orient landscape
%                 print('-dpdf', [dirFortran FigName]);
            end
        end
    else
        LeastSquaresTrixD = nan*zeros(length(UniLineNum),length(UniStationNum));
        SquidProfile(z).LeastSquaresDescent = nan*ones(size(LeastSquaresTrixD));
        % now pull out indices where they are the lowest (least squares)
        [newD, indD] = sort(LeastSquaresTrixD(:), 1, 'ascend');
        newDtrix = [newDtrix; newD];
        lsvalD(:,z) = newD(1:nn);
        lsindD(:,z) = indD(1:nn);
        lsdateD(:,z) = DatetrixD(lsindD(:,z));
    end
end
if plotOLSprofiles
    if panelnum < numberofpanels
        key = ['CalCOFI (gray) v. Squid ' num2str(TagHinke) 'Descents (pink)'];
        [h3] = suptitle(key);
        [ax,h(1)] = suplabel(['temp'],'x');
        [ax,h(2)] = suplabel('depth','y');
        set(h3,'FontSize',24, 'FontWeight', 'bold')
        set(h,'FontSize',15, 'FontWeight', 'bold')
        FigName = [num2str(TagHinke) 'olsCalCOFI_D_' num2str(b) 'hrbuff' num2str(fignum) '.pdf'];
        annotate_JS(Mfilename, gcf, FigName)
%         orient landscape % save
%         print('-dpdf', [dirFortran FigName]);
    end
end

% make lat and lon into a vector with indices instead of matrices
newLonD = [];
newLatD = [];
for k = 1:size(LontrixD,2)
    newLonD = [newLonD; LontrixD(:,k)];
    newLatD = [newLatD; LattrixD(:,k)];
end

newDateD = [];
for k = 1:size(DatetrixD,2)
    newDateD = [newDateD; DatetrixD(:,k)];
end

% plot RMS histos
figure; hold on
[n, xout] = hist(newDtrix,50);
[n2, xout2] = hist(lsvalD(:,11),xout); % for stacking
N = [n2; n; nan(1,length(n2))]'; % to make the colormap work with stacking
bar(N, 'stacked')
colormap(gray)
xlabel('RMS value ','FontSize', 14, 'FontWeight', 'bold')
ylabel('Count', 'FontSize', 14, 'FontWeight', 'bold')
title({[num2str(TagHinke) ' RMS histogram: Best' num2str(nn) 'Descents ']; 'RMS of All Squid profiles v. all CalCOFI casts'}, 'fontsize', 16, 'fontweight', 'bold');
FigName = [num2str(TagHinke) '_histogram_RMS_D_' num2str(b) 'hrbuff.pdf']; %FigName = [num2str(TagHinke) '_timeseries_distoffshore_perc.pdf'];
annotate_JS(Mfilename, gcf, FigName);
% orient landscape % save
% print('-dpdf', [dirFortran FigName]);

save -mat calcofi_tags_olsD.mat UniCoords SpatialLon SpatialLat lsindD lsvalD newLonD newLatD

%% plot both RMS histos together with ascents and descents

figure; hold on
% [n, xout] = hist([newAtrix; newDtrix],50); % STACK THIS SHIT
% [n2, xout2] = hist([lsvalA(:,11); lsvalD(:,11)],xout); % for stacking
% N = [n2; n; nan(1,length(n2))]'; % to make the colormap work with stacking
[n, xout] = hist([newAtrix],50); 
[n1, xout] = hist([newDtrix],50); 
[n2, xout2] = hist([lsvalA(:,11); lsvalD(:,11)],xout); % for stacking
N = [n2; n; n1; nan(1,length(n2))]'; % to make the colormap work with stacking
bar(N, 'stacked')
colormap(gray)
xlabel('RMS value ','FontSize', 14, 'FontWeight', 'bold')
ylabel('Count', 'FontSize', 14, 'FontWeight', 'bold')
title({[num2str(TagHinke) ' RMS histogram: Best' num2str(nn) 'Ascents and Descents ']; 'RMS of All Squid profiles v. all CalCOFI casts'}, 'fontsize', 16, 'fontweight', 'bold');
FigName = [num2str(TagHinke) '_histogram_RMS_AD_' num2str(b) 'hrbuff.pdf']; 
annotate_JS(Mfilename, gcf, FigName);
% orient landscape % save
% print('-dpdf', [dirFortran FigName]);

%% plotting!

% doesn't make sense to plot both on the same thing, since ascent and
% descent of squid are different.... this was a work in progress July 14
% calcofi_tags_olsPlots % delete this once it works in the loops above.

%% test the indices from above
% lsind is created in radiusApproach_Tags.m

if plothere
    for z = 1:length(SquidProfile)
        figure
        mapSqArea
        hold on
        set(gca, 'XLim', [-125 max(SpatialLon)], 'YLim', [min(SpatialLat)-1 37]);
        plot(UniCoords(1:end-1,1),UniCoords(1:end-1,2), '.', 'Color', [0.8 0.8 0.8]) % 1:end-1 so it's the same length as UniCasts
        plot(newLonA(lsindA(1:nn,z)),newLatA(lsindA(1:nn,z)), 'om',  'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63], 'MarkerSize', 11)
        set(gca, 'fontsize', 12, 'fontweight', 'bold')
        xlabel('Longitude','FontSize', 14, 'FontWeight', 'bold')
        ylabel('Latitude','FontSize', 14, 'FontWeight', 'bold')
        title({['5-83051 Day ' num2str(z)];['CalCOFI ' num2str(nn) ' least squares']}, 'FontSize', 18, 'FontWeight', 'bold');
        FigName = ['CalCOFI_ols' num2str(nn) '_map_' num2str(b) 'hrbuff_Day' num2str(z) '.pdf'];
        annotate_JS(Mfilename, gcf, FigName);
        %     orient landscape % save
        %     print('-dpdf', [dirCTD FigName]);
    end
end

%%

disp('Completed calcofi_tags_ols.m')
% ===== EOF [calcofi_tags_ols.m] ======
