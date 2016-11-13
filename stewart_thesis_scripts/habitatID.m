% habitatID.m
%  
% This is called by manageTagData.m 
%  
% Outside Functions Called: 
        % compare_satellite.m (J. Stewart)
        % probHabitat.m (J. Stewart)
        % mapNEPac (A. Booth)
%  
% AUTHOR    : J. Stewart jules32@gmail.com 
% DATE      : 11-Aug-2010 12:01:35  
% Revision  : 1.00  
% DEVELOPED : 7.9.0.529 (R2009b) OSX 
% FILENAME  : compareTagsSatellite.m 

Mfilename = mfilename;

%% This worked from habitatID_eachTag on 07-Mar-2011, but should not use
%% polyval because not normal!
%% Find interquartile range of each day this isn't a normal distribution...
% DO THIS AS IQR


polyOrder = 7;
X = binincr;
polyCoeftrix = [];
for i = 1:size(HinkeHistCounts,1)
    HistCountsPerc = HinkeHistCounts(i,:)/sum(HinkeHistCounts(i,:)); % before running polyval, normalize to 1. then polyCoef will be a decimal. 
    Y = HistCountsPerc;
    [polyCoef,polystats] = polyfit(X,Y,polyOrder); % DOES NOT work with NaNs!!!!!!!!!!  
    polyCoeftrix = [polyCoeftrix; polyCoef];
end



%% Code from Dave Foley (11-Aug-2010) from methods used in Hinke et al.   
% exact commands from D. Foley (11-Aug-2010) saved in habitatID.m

%% This runs polyfit, although this isn't a normal distribution...

polyOrder = 7;
X = binincr;
polyCoeftrix = [];
for i = 1:size(HinkeHistCounts,1)
    HistCountsPerc = HinkeHistCounts(i,:)/sum(HinkeHistCounts(i,:)); % before running polyval, normalize to 1. then polyCoef will be a decimal. 
    Y = HistCountsPerc;
    [polyCoef,polystats] = polyfit(X,Y,polyOrder); % DOES NOT work with NaNs!!!!!!!!!!  
    polyCoeftrix = [polyCoeftrix; polyCoef];
end

% % plot
figure
h2 = bar(X, Y,'BarWidth',1);
set(h2,'facecolor',[0.7529    0.7529    0.7529],'EdgeColor',[0.5020    0.5020    0.5020],'LineWidth',1.5)
hold on
[polyCoef,polystats] = polyfit(X,Y,polyOrder); % DOES NOT work with NaNs!!!!!!!!!! % from A. Booth's Example_LinearFit.m. See help polyfit for a good discussion of outputs
slope = polyCoef(1);
yinter = polyCoef(2);

Xvec = floor(min(binincr)):0.1:floor(max(binincr));
polyFit = polyval(polyCoef,Xvec);
% stats = regstats(Y,X,'poly',{'rsquare'});
% Rsqr = stats.rsquare;
plot(X,Y,'*')
hold on
plot(Xvec,polyFit,'r:', 'LineWidth', 4)
hold off
xlabel('Temperature','FontSize', 14, 'FontWeight', 'bold')
ylabel('Percent','FontSize', 14, 'FontWeight', 'bold')
title({['Combined Tag Temperature Histograms: Shallow (<' num2str(ShallowCutoff) 'm),']; [num2str(polyOrder) 'th Order Polynomial']}, ...
    'FontSize', 14, 'FontWeight', 'bold');  
FigName = [num2str(TagHinke) 'TempHistosShallow' num2str(ShallowCutoff) 'mPolyFit_' num2str(polyOrder) '.pdf'];
annotate_JS(Mfilename, gcf, FigName)
% orient landscape % save
% print('-dpdf', [dirHinke FigName]);

%% Code from Dave Foley (11-Aug-2010) from methods used in Hinke et al.   
% exact commands from D. Foley (11-Aug-2010) saved in habitatID.m

for i = 1:length(stime)
    figure
    sstd = squeeze(sst_3D(i,:,:));  
    pcolor(slon,slat,sstd)
    ind = sstd < 13 | sstd > 15; % turn anything outside of the temp boundaries in sstd into nan
    sstd(ind) = nan;
    phui = polyval(polyCoeftrix(i,:),sstd); % probable habitat utilization index. % polyval can also do error estimates
    ind = phui < 0;
    phui(ind) = 0; % change any negative phui values to 0 [there shouldn't be any negatives at this point, maybe don't need this line]
    ind = isnan(phui);
    phui(ind) = 0; % turn any nan phui values to 0
    phui = phui./max(max(phui)); % normalize so scalebar will be 0 to 1
    pcolor(slon,slat,phui) % draw the figure
    Hcbar = colorbar('location', 'eastoutside');
    colormap(gray) % colormap should be from 0 to 1.
    ylabel(Hcbar,'P(Habitat: Suitable Surface Temperature)  ', 'FontSize',12, 'FontWeight', 'bold');
    title([num2str(TagHinke) ' ' datestr(stime(i), 1) ], 'FontSize',20, 'FontWeight', 'bold')
    % pause
    hold on % now plot the coastline
    mapNEPac(Locations, NEPacBuffer, []) %1.5 works when LocationsCA is 3 Monterey 2009 tags 
    hold on
    plot(LocTagCA_Hinke(:,1),LocTagCA_Hinke(:,2),'r.', 'markersize', 20)
    FigName = [num2str(TagHinke) '_HabitatID_' num2str(ShallowCutoff) 'm_' datestr(stime(i), 'yyyy_mmm_dd') '.pdf'];
    annotate_JS(Mfilename, gcf, FigName)
%     orient landscape % save
%     print('-dpdf', [dirHinke FigName]);
end
% %% trying to work with IQR and wind. Work in progress. 
% 
% QuantileChoice = [2 6]; %[1 7]; % in tagDailyHistos_Shallow,  QuantileSet = [0.025, 0.05, 0.25 0.5, 0.75, 0.95, 0.975];
% 
% for i = 1:size(TempShallowQuantile,1)
%     
%     HistCountsPerc = HinkeHistCounts(i,:)/sum(HinkeHistCounts(i,:)); % before running polyval, normalize to 1. then polyCoef will be a decimal. 
%     Y = HistCountsPerc;
%     [polyCoef,polystats] = polyfit(X,Y,polyOrder); % DOES NOT work with NaNs!!!!!!!!!!  
%     polyCoeftrix = [polyCoeftrix; polyCoef];
% end
% 
% for i = 1:length(stime)
%     figure
%     sstd = squeeze(sst_3D(i,:,:));  
%     pcolor(slon,slat,sstd)
%     ind = sstd < TempShallowQuantile(i,QuantileChoice(1)) | sstd > TempShallowQuantile(i,QuantileChoice(2)); 
%     sstd(ind) = 0; % setup to turn anything outside of the temp boundaries in sstd into white
%     phui = polyval(TempShallowQuantile(i,QuantileChoice),sstd); % probable habitat utilization index. % polyval can also do error estimates
%     ind = phui < 0;
%     phui(ind) = 0; % change any negative phui values to 0 [there shouldn't be any negatives at this point, maybe don't need this line]
%     ind = isnan(phui);
%     phui(ind) = 0; % turn any nan phui values to 0. Must change this to grey for clouds/irretreivable data
%     phui = phui./max(max(phui)); % normalize so scalebar will be 0 to 1
%     pcolor(slon,slat,phui) % draw the figure
%     Hcbar = colorbar('location', 'eastoutside');
%     colormap(gray) % colormap should be from 0 to 1.
%     ylabel(Hcbar,'P(Habitat: Suitable Surface Temperature)  ', 'FontSize',12, 'FontWeight', 'bold');
%     title([num2str(TagHinke) ' ' datestr(stime(i), 1) ], 'FontSize',20, 'FontWeight', 'bold')
%     % pause
%     hold on % now plot the coastline
%     mapNEPac(Locations, NEPacBuffer, []) %1.5 works when LocationsCA is 3 Monterey 2009 tags 
%     hold on
%     plot(LocTagCA_Hinke(:,1),LocTagCA_Hinke(:,2),'r.', 'markersize', 20)
%     FigName = [num2str(TagHinke) '_HabitatID_' num2str(ShallowCutoff) 'm_' datestr(stime(i), 'yyyy_mmm_dd') '.pdf'];
%     annotate_JS(Mfilename, gcf, FigName)
% %     orient landscape % save
% %     print('-dpdf', [dirHinke FigName]);
% end
% % 
% % 
% % %% wind analysis
% % 
% % Windtrix = [];
% % for i = 1:length(wtimeZ)
% %     figure
% %     windd = squeeze(wind_3D(i,:,:)); % this is just for the third day.   
% % %     reshape with the dimensions of the new matrix
% %     pcolor(wlonZ,wlatZ,windd)
% %     ind = windd < 13 | windd > 15; % turn anything outside of the temp boundaries in sstd into nan
% %     windd(ind) = nan;
% %     Windtrix = [Windtrix windd];
% %     
% % %     phui = polyval(polyCoeftrix(i,:),windd); % probable habitat utilization index. % polyval can also do error estimates
% % %     ind = phui < 0;
% % %     phui(ind) = 0; % change any negative phui values to 0 [there shouldn't be any negatives at this point, maybe don't need this line]
% % %     ind = isnan(phui);
% % %     phui(ind) = 0; % turn any nan phui values to 0
% % %     phui = phui./max(max(phui)); % normalize so scalebar will be 0 to 1
% % %     pcolor(wlonZ,wlatZ,phui) % draw the figure
% % %     Hcbar = colorbar('location', 'eastoutside');
% % %     colormap(gray) % colormap should be from 0 to 1.
% % %     ylabel(Hcbar,'P(Habitat: Suitable Surface Temperature)  ', 'FontSize',12, 'FontWeight', 'bold');
% % %     title([num2str(TagHinke) ' ' datestr(stime(i), 1) ], 'FontSize',20, 'FontWeight', 'bold')
% % %     % pause
% % %     hold on % now plot the coastline
% % %     mapNEPac(Locations, NEPacBuffer, []) %1.5 works when LocationsCA is 3 Monterey 2009 tags 
% % %     hold on
% % %     plot(LocTagCA_Hinke(:,1),LocTagCA_Hinke(:,2),'r.', 'markersize', 20)
% % %     FigName = [num2str(TagHinke) '_HabitatID_' num2str(ShallowCutoff) 'm_' datestr(stime(i), 'yyyy_mmm_dd') '.pdf'];
% % %     annotate_JS(Mfilename, gcf, FigName)
% % %     orient landscape % save
% % %     print('-dpdf', [dirHinke FigName]);
% % end


%%

% This data is not good for fitting a polynomial. So document this, but
% we'll probably do something just based on bins. So make the bins coarser.
% With binsize of .5 could probably say that it was normally distributed. 

% run a median filter %%%TODO******
% normalize TempHist before fitting a curve

TempHistPerc = TempHist/sum(TempHist); % before running polyval, normalize to 1. then polyCoef will be a decimal. 

X = binT;
Y = TempHistPerc;
polyOrder = 7;

% plot
figure
h2 = bar(X, Y,'BarWidth',1);
set(h2,'facecolor',[0.7529    0.7529    0.7529],'EdgeColor',[0.5020    0.5020    0.5020],'LineWidth',1.5)
hold on

% from A. Booth's Example_LinearFit.m. See help polyfit for agood disssion of outputs
[polyCoef,polystats] = polyfit(X,Y,polyOrder); % DOES NOT work with NaNs!!!!!!!!!!
slope = polyCoef(1);
yinter = polyCoef(2);

Xvec = floor(min(TempShallow)):0.1:floor(max(TempShallow));
polyFit = polyval(polyCoef,Xvec);
% stats = regstats(Y,X,'poly',{'rsquare'});
% Rsqr = stats.rsquare;
plot(X,Y,'*')
hold on
plot(Xvec,polyFit,'r:', 'LineWidth', 4)
hold off
xlabel('Temperature','FontSize', 14, 'FontWeight', 'bold')
ylabel('Percent','FontSize', 14, 'FontWeight', 'bold')
title({['Combined Tag Temperature Histograms: Shallow (<' num2str(ShallowCutoff) 'm),']; [num2str(polyOrder) 'th Order Polynomial']}, ...
    'FontSize', 14, 'FontWeight', 'bold');  
FigName = ['TempHistosCombinedShallow' num2str(ShallowCutoff) 'mPolyFit_' num2str(polyOrder) '.pdf'];
annotate_JS(Mfilename, gcf, FigName)
% orient landscape % save
% print('-dpdf', FigName);

% do counts for histograms, not percents. percents will normalize so that it's not for a particular animal so you don't give the mx one more weight
% ultimately, look at day+night together, but look at them separately too. 
% decrease bin size so that can use higher order polynomial (see error
% message)
% look at COASTWATCH browser, and put the temperature range in the
% narrowest band of this (eg 12-14°C) and it will show you where it's too
% hot, where it's too cold in terms of the probability of finding squid. 


%% Read in SST data for variables defined in get_stewart
% minlon = -118;
% maxlon = -123.5;
% minlat = 32;
% maxlat = 37;
% 
% tmin = datenum(2009, 9, 30); % deployment of tag 83046_09
% tmax = datenum(2009, 11, 24); % popoff of tag 83052_09


load('get_stewartOutput.mat') % created in get_stewart: uses xtractomatic

%% Tag Information inputed in CompareTags.m

Locations = [LocTagDeployCA; LocTagPopUpCA];
Locations = [Locations(:,2) Locations(:,1)];

%     Labels = [
%        % 1;
%         2;
%         3;
%         4;
%         5;
%       %  1.1;
%         2.1;
%         3.1;
%         4.1;
%         5.1
%         ];

%% Code from Dave Foley (11-Aug-2010) from methods used in Hinke et al.   

for i = 1:2%length(stime)
    figure
    sstd = squeeze(sst_3D(i,:,:)); % this is just for the third day.   
    pcolor(slon,slat,sstd)
    ind = sstd < 13 | sstd > 15; % turn anything outside of the temp boundaries in sstd into nan
    sstd(ind) = nan;
    phui = polyval(polyCoef,sstd); % probable habitat utilization index. % polyval can also do error estimates
    ind = phui < 0;
    phui(ind) = 0; % change any negative phui values to 0 [there shouldn't be any negatives at this point, maybe don't need this line]
    ind = isnan(phui);
    phui(ind) = 0; % turn any nan phui values to 0
    phui = phui./max(max(phui)); % normalize so scalebar will be 0 to 1
    pcolor(slon,slat,phui) % draw the figure
    Hcbar = colorbar('location', 'eastoutside');
    colormap(gray) % colormap should be from 0 to 1.
    ylabel(Hcbar,'P(Habitat: Suitable Surface Temperature)  ', 'FontSize',12, 'FontWeight', 'bold');
    title(datestr(stime(i), 1), 'FontSize',20, 'FontWeight', 'bold')
    % pause
    hold on % now plot the coastline
    mapNEPac(Locations,0.2, 1:8) %no labels % FILL
    FigName = ['HabitatID_' datestr(stime(i), 'yyyy_mmm_dd') '.pdf'];
%    annotate_JS(Mfilename, gcf, FigName) % issues
%     orient landscape % save
%     print('-dpdf', FigName);
end




% % this works
% figure
% sstd = squeeze(sst_3D(3,:,:)); % this is just for the third day. 
% pcolor(slon,slat,sstd)
% 
% phui = polyval(polyCoef,sstd); % probable habitat utilization index
% %pcolor(slon,slat,phui);
% ind = phui < 0;
% phui(ind) = 0;
% ind = sstd < 13 | sstd > 15;
% sstd(ind)=nan;
% ind = isnan(phui);
% phui(ind) = 0;
% phui = phui./max(max(phui)); % normalize so scalebar will be 0 to 1
% pcolor(slon,slat,phui)
% colorbar
% colormap(gray)

%% exact commands from Foley
% load get_stewartOutput.mat
% who
% whos
% stime(1)
% datestr(stime(1))
% datestr(stime(2))
% datestr(stime(3))
% sstd = sst_3D(3,:,:);
% pcolor(slon,slat,sstd)
% sstd = squeeze(sst_3D(3,:,:));
% pcolor(slon,slat,sstd)
% w
% help polyfit
% help polyval
% whos
% phui = polyval(polyFit,sstd);
% ind = find(phui<0);
% phui(ind)=0;
% figure
% pcolor(slon,slat,phui)
% colorbar
% min(phui)
% ls
% ind = find(sstd<11 | sstd>15);
% sstd(ind)=nan;
% phui=polyval(polyFit,sstd);
% pcolor(slon,slat,phui)
% colorbar
% pcolor(slon,slat,sstd)
% colormap
% colorbar
% polyval(polyFit,13)
% polyFit
% polyval(polyCoef,13)
% phui=polyval(polyCoef,sstd)
% pcolor(slon,slat,phui)
% ind=find(isnan(phui));
% phui(ind)=0
% pcolor(slon,slat,phui)
% colorbar
% colormap(gray)
%% 
  
disp('Completed compareTagsSatellite.m') 
% ===== EOF [compareTagsSatellite.m] ======  
