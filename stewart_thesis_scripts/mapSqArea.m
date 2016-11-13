    % mapSqArea.m
%  this is called from c_bathymetry.m by J. Stewart. This is modified from
%  mapNEPac.m, but that wasn't working properly so this is isolated. 
% Outputs: 
%  
% Outside Functions Called: 
%  
% AUTHOR    : J. Stewart jules32@gmail.com 
% DATE      : 21-Jan-2011 16:23:39  
% Revision  : 1.00  
% DEVELOPED : 7.9.0.529 (R2009b) OSX 
% FILENAME  : mapSqArea.m 

  %% get extent of points
Buffer = 10; % this determines the extent of the map, will be Buffer*100% bigger than range of points
LatExtrema = [nanmin(Locations(:,1)), nanmax(Locations(:,1))];
LonExtrema = [nanmin(Locations(:,2)), nanmax(Locations(:,2))];
LatBuffer = range(LatExtrema) * Buffer; % this determines the extent of the map
LonBuffer = range(LonExtrema) * Buffer; % will be 20% bigger than range of points
if LatBuffer == 0 %if only 1 location point
    LatBuffer = 0.5 * Buffer; % this determins the extent of the map
    LonBuffer = 0.5 * Buffer; % will be 20% bigger than range of points
end
latrange = [LatExtrema(1)-LatBuffer, LatExtrema(2)+LatBuffer]; % map range
lonrange = [LonExtrema(1)-LonBuffer, LonExtrema(2)+LonBuffer];

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

%% Plot
hold on
set(gcf,'Color','w')
plot(ngdc_NEPac(:,1),ngdc_NEPac(:,2),'Color',[0.5 0.5 0.5]);
set(gca,'DataAspectRatio',[111 cos(37*pi/180)*111 1],'YLim', latrange, 'XLim', lonrange);
hold on
plot(Locations(:,2),Locations(:,1),'w.')
% scatter(Locations(:,2),Locations(:,1),40,1:length(Locations(:,1)),'filled') %for numerical ids
hold off

%% 
  
disp('Completed mapSqArea.m') 
% ===== EOF [mapSqArea.m] ======  
