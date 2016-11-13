% mappingPostFortran_dist.m
%  
% Called From: habitatID_eachTag.m % J. Stewart
  
% Description: calculates minimum distances from the coast, the shelf, and
% also distances between points. Makes figures. 
%  
% Outside Functions Called: 
% load get_stewart_SSTbathy_83051.mat created in get_stewart.m % J. Stewart
% distance_haversine % J. Stewart
%  
% AUTHOR    : J. Stewart jules32@gmail.com 
% DATE      : 18-Jul-2011 13:14:44  
% Revision  : 1.00  
% DEVELOPED : 7.9.0.529 (R2009b) OSX 
% FILENAME  : mappingPostFortran_dist.m 

Mfilename = mfilename; 

%%

% m/s calcs:
swimmax*...(km/d)
1000/...(m/km)
24/...(d/hr)
60/...(hr/min)
60 %(min/sec)

% bodylengths calcs: to get bodylengths/sec
swimmax*...(km/d)
1000/...(m/km)
0.79/...(bodylengths/m)
24/...(d/hr)
60/...(hr/min)
60 %(min/sec) 

taggingDistHaversine = distance_haversine(LocTagCA_Hinke(1,1), LocTagCA_Hinke(1,2), LocTagCA_Hinke(2,1),LocTagCA_Hinke(2,2))  

%%
Locations = [34.4,-122;34.5,-122.1];
% figure

%% get extent of points

latrange = [31, 38];
lonrange = [-126, -116];

%% crop ngdc_NEPac

disp('Looking for ngdc_NEPac.dat from http://rimmer.ngdc.noaa.gov: Range: 46.3->22.5N and -125->-107W,  matlab format, World Data Bank II')
file = 'ngdc_NEPac.dat'; %file from http://rimmer.ngdc.noaa.gov/matlab format "World Data Bank II (designed for 1:2,000,000)"
if exist(file,'file')
    disp('ngdc_NEPac.dat file found')
	load ngdc_NEPac.dat 
else
    file = uigetfile('*.dat','Please locate ngdc_NEPac.dat or create file: <http://rimmer.ngdc.noaa.gov/> 46.3->22.5N and -125->-107W, "World Data Bank II (designed for 1:2,000,000)", matlab format).');
    load(file)
end    

% crop ngdc_NEPac
ngdc_NEPac_uncrop = ngdc_NEPac;

for k = 1:length(ngdc_NEPac_uncrop)% change anything outside of map range to NaN
    if (ngdc_NEPac_uncrop(k,1)>lonrange(2)) || (ngdc_NEPac_uncrop(k,1)<lonrange(1));
            ngdc_NEPac_uncrop(k,:)=NaN;
    elseif (ngdc_NEPac_uncrop(k,2)>latrange(2)) || (ngdc_NEPac_uncrop(k,2)<latrange(1));
            ngdc_NEPac_uncrop(k,:)=NaN;
    end
end
clear k
crop = find(isfinite(ngdc_NEPac_uncrop(:,1))); % row index of data points
ngdc_NEPac = ngdc_NEPac_uncrop(crop,:); % remove all rows with NAN

%% Prevent island effect (drawing lines between islands and mainland)
%     hist(abs(difflon),[-.5:0.005:3]) %to see range in distances to choose threshhold 
threshhold = 0.04; %if distance between point is greater than threshold change number to NaN
difflon = diff(ngdc_NEPac(:,1));     
ngdc_NEPac(abs(difflon)>threshhold,:) = NaN;
difflat = diff(ngdc_NEPac(:,2));
ngdc_NEPac(abs(difflat)>threshhold,:) = NaN;

% remove inland islands % J. Stewart 18-Jul-11
ngdc_NEPac(ngdc_NEPac(:,1)>-121.2 &ngdc_NEPac(:,2)>36) = NaN; % remove lakes
ngdc_NEPac(ngdc_NEPac(:,1)<-119.35 &ngdc_NEPac(:,2)<34.2) = NaN; % remove channel islands
ngdc_NEPac(ngdc_NEPac(:,1)<-118.2 &ngdc_NEPac(:,2)<33.6) = NaN;% remove channel islands

%% Plot

% hold on
% set(gcf,'Color','w')
% plot(ngdc_NEPac(:,1),ngdc_NEPac(:,2),'Color',[0.5 0.5 0.5]);
% set(gca,'DataAspectRatio',[111 cos(37*pi/180)*111 1],'YLim', latrange, 'XLim', lonrange);

%% setup

countProb1 = []; % count for how many pixels with probability==1
for d = 1:length(probDayUni)
    dlog = probDayUni(d) == probDay; % called 'good' in R code
    plog = probProb(dlog) > 0; % called 'good' in R code
    countProb1 = [countProb1; sum(plog)]; % count
end

%% calculate the range of minimums. So find the minimum distance from the coast ("offshore")
% point to the coast and then for each day have a range of those minimum distances.  
  
% OffshoreTrix = nan(max(countProb1),length(probDayUni)); % set up for min, mean, max distance for each day
% wmean = [];
% figure; hold on
% for d = 1:length(probDayUni)
%     dlog = probDayUni(d) == probDay;
%     plog = probProb > 0;
%     dp = dlog+plog;
%     dplog = dp == 2;
%     x = probLon(dplog);
%     y = probLat(dplog);
%     
%     disttrix = [];
%     for i = 1:length(x)
%         [dist_km] = distance_haversine(x(i), y(i), ngdc_NEPac(:,1), ngdc_NEPac(:,2));
%         disttrix = [disttrix; min(dist_km)];
%     end
%     OffshoreTrix(1:length(x),d) = disttrix;
%     wmean = [wmean; sum(disttrix.*probProb(dplog))];
%     scatter(repmat(d-1,length(x),1),disttrix,probProb(dplog)*1000,'markeredgecolor', 'k')
% end
% 
% disp('mean daily dist traveled offshore')
% mean(diff(wmean))
% 
% plot(0:length(probDayUni)-1,wmean, 'k', 'linewidth', 1.5)
% xlabel('Day ','FontSize', 14, 'FontWeight', 'bold')
% ylabel({'Distance offshore (km), weighted by probability (size) '; 'weighted mean (black line)'}, 'FontSize', 14, 'FontWeight', 'bold')
% title([num2str(TagHinke) ' distance offshore '], 'fontsize', 16, 'fontweight', 'bold');
% set(gca, 'XLim', [-1 Deployment_day], 'fontsize', 12, 'fontweight', 'bold')
% % h1 = gca;
% % h2 = axes('Position', get(h1, 'Position'));
% % plot(0:length(probDayUni)-1, countProb1/sum(msk)*100, '-', 'color', [0.8 0.8 .8], 'linewidth', 2) % put this on other axis
% % set(h2, 'YAxisLocation', 'right', 'XTickLabel', [], 'Color', 'none', 'fontsize', 12, 'fontweight', 'bold', 'XGrid','off','YGrid','off','Box','off');
% % set(h2, 'XLim', get(h1, 'XLim'), 'YLim', [0 25], 'Layer', 'top');
% % ylabel('percent of non-zero probability tagging area', 'fontsize', 16, 'fontweight', 'bold', 'color', [0.8 0.8 .8]);
% % legend('big', 'bigger')
% FigName = [num2str(TagHinke) '_timeseries_distoffshore.pdf']; %FigName = [num2str(TagHinke) '_timeseries_distoffshore_perc.pdf'];
% annotate_JS(Mfilename, gcf, FigName);
% orient landscape % save
% % print('-dpdf', [dirFortran FigName]);

%% distance between days

% from Ricardo:
% I shouldn't be responding at this hour, I am no longer a grad student :-) But since you are going at it, who am I to stop you!
% I think the maximum distance that you are computing also has the problem of trying to connect points that cannot be connected, 
% because they are too far apart - you're getting distances well beyond 55km. Hmm. Here's a late night alternative:
% 
% For each point P with probability>0 on day d, you calculate the distances to each point with probability>0 on day d+1, 
% and exclude those that exceed the critical value (55km). You perform a weighted average of the acceptable ones, 
% based on the probabilities at day d+1 (don't forget the normalizing denominator*), and so you get the average displacement of the squid, 
% supposing she was at that point P on day d. You do this operation for all the points on day d, meaning that it shouldn't be too much information to plot. 
% Finally, you do the weighted average of those values, by using the probabilities from day d (this average does not require normalization, 
% because the sum of the probabilities for any given day equals 1).
% 
% * the normalizing denominator is the sum of the probabilities of the cells at day d+1 that are within the range of point P at day d. In other words:
% Weighted average = [ Sum (for all points within range) of distance * probability at d+1 ] / Sum (for all points within range) of probability at d+1

DailyTrix_sw = zeros(max(countProb1),length(probDayUni)); % sw means swimmax
wmean_sw = [];
figure; hold on
for d = 1:length(probDayUni)
    if d > 1
        % id today
        dlog = probDayUni(d-1) == probDay;
        plog = probProb > 0;
        dp = dlog+plog;
        dplog = dp == 2;
        prob = probProb(dplog);
        x = probLon(dplog);
        y = probLat(dplog);
        
        % id tomorrow
        dlog1 = probDayUni(d) == probDay;
        dp1 = dlog1+plog;
        dplog1 = dp1 == 2;
        prob1 = probProb(dplog1);
        x1 = probLon(dplog1);
        y1 = probLat(dplog1);
        
        disttrix_sw = [];
        for i = 1:length(x1)
            [dist_km1] = distance_haversine(x1(i), y1(i), x, y);
%             [dist_deg] = distance(x1(i), y1(i), x, y); % distance doesn't get it where it needs to be: only 545 km total dist
%             dist_km1 = deg2km(dist_deg);
            swind = find(dist_km1<swimmax);
            wm = nansum(dist_km1(swind).*prob(swind))/nansum(prob(swind));
            disttrix_sw = [disttrix_sw; wm];
        end
        DailyTrix_sw(1:length(x1),d) = disttrix_sw;
        wmean_sw = [wmean_sw; nansum(probProb(dplog1).*disttrix_sw)];
        scatter(repmat(d-1,length(x1),1),disttrix_sw,probProb(dplog1)*1000,'markeredgecolor', 'k')
%     else
%         scatter(d-1,0,1000,'markeredgecolor', 'k')
    end
end

disp('mean velocity (from daily weighted means) ±SD:')
disp(mean(wmean_sw)) % mean weighted traveling speed
disp(std(wmean_sw))
disp('max velocity (from daily weighted means):')
disp(max(wmean_sw)) % max daily traveling speed
disp('total distance traveled (from daily weighted means):')
disp(sum(wmean_sw)) % mean weighted traveling speed

% plot(0:length(probDayUni)-1,[0; wmean_sw], 'k', 'linewidth', 1.5)
plot(1:length(probDayUni)-1,wmean_sw, 'k', 'linewidth', 1.5)
xlabel('Day ','FontSize', 14, 'FontWeight', 'bold')
ylabel({'Distance from previous day (km), weighted by probability (size) '; 'weighted mean (dark line)'}, 'FontSize', 14, 'FontWeight', 'bold')
title({[num2str(TagHinke) ' distance between each day ']; ...
    %     'and percent of non-zero probability values in tagging area (ocean only) ' ...
    }, 'fontsize', 16, 'fontweight', 'bold');
set(gca, 'XLim', [-1 ceil(Deployment_day)], 'fontsize', 12, 'fontweight', 'bold')%, 'XTickLabel', [])
h1 = gca;
% h2 = axes('Position', get(h1, 'Position'));
% plot(0:length(probDayUni)-1, countProb1/sum(msk)*100, '-', 'color', [0.8 0.8 .8], 'linewidth', 2) % put this on other axis
% set(h2, 'YAxisLocation', 'right', 'XTickLabel', [], 'Color', 'none', 'fontsize', 12, 'fontweight', 'bold', 'XGrid','off','YGrid','off','Box','off');
% set(h2, 'XLim', get(h1, 'XLim'), 'YLim', [0 25], 'Layer', 'top');
% ylabel('percent of non-zero probability tagging area', 'fontsize', 16, 'fontweight', 'bold', 'color', [0.8 0.8 .8]);
FigName = [num2str(TagHinke) '_timeseries_distdaily.pdf']; % FigName = [num2str(TagHinke) '_timeseries_distdaily_perc.pdf'];
annotate_JS(Mfilename, gcf, FigName);
% orient landscape % save
% print('-dpdf', [dirFortran FigName]);

%% calc dist to shelf contour, not coast

% load -mat CenCalBathy_5706XIYI.mat %Tried this way, from c_bathymetry. Instead, do it from get_stewart.
% or try this way
% load -mat  get_stewart_SSTbathy_83051.mat %created in get_stewart.m % J. Stewart
 ShelfBreak = -500; % set contour. map can help you ID which ones. In order to not get a glitch with mini islands, do the deep part of the slope.
% Interval = 0.001; % make the interval super small since I can't figure out how to just map one contour line
% 
% figure
% hold on
% [C, h] = contour(BLON-360,BLAT,-bathy); % so can see some of the contour lines
% clabel(C,h);
% plot(LocTagCA_Hinke(:,1), LocTagCA_Hinke(:,2), 'sg', 'markersize', 6)
% 
% figure
% hold on
% [C, h] = contour(BLON-360,BLAT,-bathy, ShelfBreak-Interval:Interval:ShelfBreak); % just at specified shelf break
% clabel(C,h);
% plot(LocTagCA_Hinke(:,1), LocTagCA_Hinke(:,2), 'sg', 'markersize', 6)
% C(:,C(1,:)<=ShelfBreak) = [];
% XX = C(1,2:end);
% YY = C(2,2:end);
% 
% 
% threshhold = 300; %if distance between point is greater than threshold change number to NaN
% difflon = diff(XX); 
% difflat = diff(YY);
% difflonind = find(abs(difflon)>1);
% difflatind = find(abs(difflat)>10);
% 
% XX(difflonind-1) = NaN;
% YY(difflatind-1) = NaN;
% 
% figure
% hold on
% plot(XX,YY, '.')
% plot(LocTagCA_Hinke(:,1), LocTagCA_Hinke(:,2), 'sg', 'markersize', 6)
% set(gca,'XLim', [-124 -116], 'YLim', [32 37])

% Do it this way with a text file Patrick cleaned up in ArcGIS. there were
% too many little islands for the 500 m contour. 

filein = '/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/etopo2_500m.txt';
fid = fopen(filein,'rt');
temp = textscan(fid, '%f %f', 'delimiter','\t', 'headerlines', 1);
fclose(fid);
[XX YY] = temp{:}; % these are Lon and Lat for contour500


ContourTrix = nan(max(countProb1),length(probDayUni)); % set up for min, mean, max distance for each day
wmeanC = []; % Contour
figure; hold on
for d = 1:length(probDayUni)
    dlog = probDayUni(d) == probDay;
    plog = probProb > 0;
    dp = dlog+plog;
    dplog = dp == 2;
    x = probLon(dplog);
    y = probLat(dplog);
    
    disttrixC = [];
    onshelfC = [];
    for i = 1:length(x)
        %[dist_kmC, az] = distance(x(i), y(i), XX, YY); % do not use. would have to use deg2km too, but gives funky results.
        [dist_kmC] = distance_haversine(x(i), y(i), XX, YY);
        onshelf = (x(i) - XX) > 0;
        XXmindist = XX(find(dist_kmC==min(dist_kmC)));
        YYmindist = YY(find(dist_kmC==min(dist_kmC)));
        if x(i)-XXmindist(1) > 0 | y(i)-YYmindist(1) > 0                
            disttrixC = [disttrixC; (min(dist_kmC))];
        else
            disttrixC = [disttrixC; (min(dist_kmC)*-1)];
        end
    end
    if d == length(probDayUni) && TagHinke == 83051% last day of Tag 5 (83051) is off the shelf in SoCal Bight
        disttrixC = disttrixC*-1;
    end
    ContourTrix(1:length(x),d) = disttrixC;
    wmeanC = [wmeanC; sum(disttrixC.*probProb(dplog))];
    scatter(repmat(d-1,length(x),1),disttrixC,probProb(dplog)*1000,'markeredgecolor', 'k')
end
plot(0:length(probDayUni)-1,wmeanC, 'k', 'linewidth', 1.5)
xlabel('Day ','FontSize', 14, 'FontWeight', 'bold')
ylabel({['Distance (km) from ' num2str(ShelfBreak*-1) ' m contour, weighted by probability (size) ']; 'weighted mean (black line)'}, 'FontSize', 14, 'FontWeight', 'bold')
title([num2str(TagHinke) ' distance from shelf (' num2str(ShelfBreak*-1) ' m contour)  '], 'fontsize', 16, 'fontweight', 'bold');
set(gca, 'XLim', [-1 Deployment_day], 'fontsize', 12, 'fontweight', 'bold')
FigName = [num2str(TagHinke) '_timeseries_distContour' num2str(ShelfBreak*-1) '.pdf']; %FigName = [num2str(TagHinke) '_timeseries_distoffshore_perc.pdf'];
annotate_JS(Mfilename, gcf, FigName);
orient landscape % save
% print('-dpdf', [dirFortran FigName]);

disp('direct dist from shelf (day 1 and end): ')
disp(([wmeanC(1) wmeanC(end)])) % max weighted dist from shelf (it's min because these are negative values
disp('max dist from shelf (from daily weighted means): ')
disp(min(wmeanC)) % max weighted dist from shelf (it's min because these are negative values

% for excel sheet
% 

%% 
  
disp('Completed mappingPostFortran_dist.m') 
% ===== EOF [mappingPostFortran_dist.m] ======  
