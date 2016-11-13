% get_stewart.m  

%this is called by manageTagData

% it will access SST data and also wind data. You need to enter the days
% that you prefer. 

% program to read marlin data (full data set) and xtract enviro data
% follow get_musyl_all.m (D Foley) as a template 

% Tag data I have: 
% 1 tag from October 2008, Cordell Bank. (time series)
% 3 tags from October-November 2009, Monterey (2 time series, 1 PDT)
% 1 tag from November 2009, Monterey-Mexico (time series)


% use xtracto_3D_bdap to get a 3 dimentional block looking at SST across
% space and time. Do this in 2 different analyses: 1 for 2008 and 1 for
% 2009. **May have to change this into 3 different anaylses (Cordell,
% Monterey, Monterey-Mexico). 

% modified 08-Mar-2011 to incorporate bathymetry: from D. Foley's
% make_bathy_mask.m

cd '/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/_CA_Tags/'
% get rid of get_stewart_bathy.m

%% setup area to download, manually assign coordinates and dates

% this is based on 4 Monterey tags (2009 analysis)
xmin = 360-126;
xmax = 360-118; % 118
ymin = 31; %31 amd 38; 
ymax = 38; % switch it to 40 for 2008 tag

%¡¡¡¡ if you change the dates to deal with a different tag, change the filename below too. 
tmin = datenum(2009, 11, 05); % GMT deployment date
tmax = datenum(2009, 11, 06); % popoff date +1 or 2. check to see that stime matches up with the dates you want for the tag. 
tpos = [tmin tmax];
%% SST

% call to get SST: use xtracto_3D_bdap to get the time and space and SST you want
[sst slon slat stime] = xtracto_3D_bdap([xmin xmax],[ymin ymax],tpos,'TBAssta5day');%3°C
datestr(stime)

% assemble call to erddap for etopo2 (sandwell and smith)

% sst_3D = squeeze(sst); % this gives me [time lat lon]
% [nt ny nx] = size(sst_3D);

% animation. could try this for temperature across space with depth increments. Compare depth temps across swaths. (B. Wells)
% figure
% for i =1:nt,
%     pcolor(slon,slat,squeeze(sst_3D(i,:,:)));
%     % label with date!!!!
%     pause
% end

% save -mat get_stewart_SST_wholearea.mat %slat slon sst sst_3D stime %nt
% nx ny


%% BATHYMETRY
% assemble call to erddap for etopo2 (sandwell and smith)

% default URL for NMFS/SWFSC/ERD  THREDDS server
urlbase='http://coastwatch.pfel.noaa.gov/coastwatch/CWBrowserWW360.jsp';

% text string for data retrieval call
bobcall = strcat(urlbase,'?get=bathymetryData',...
    '&minlon=',num2str(xmin),'&maxlon=',num2str(xmax),...
    '&minlat=',num2str(ymin),'&maxlat=',num2str(ymax),...
    '&filetype=.mat')

% extract data array and import to Matlab depending on structure - a mini xtracto
varname = 'BAthym';
fileout='tmp.mat';
urlwrite(bobcall,fileout);
load('-MAT',fileout);
eval(strcat('etopo2=',varname,';'));
blon =lon;
blat=lat;

%%
% make into a 2-Degree grid and correct for sign (depth v. altitude)
bathy = -1*squeeze(etopo2);

% gonna interpolate, so set all negative value (land) to zero
ind = find(bathy<0);
bathy(ind)=0;

% define SST grid
[SLON SLAT] = meshgrid(slon,slat);

% set up bathy array on SST grid
[BLON BLAT] = meshgrid(blon,blat); 
bathy2D = griddata(BLON(:),BLAT(:),bathy(:),SLON,SLAT);


%% SAVE everything

% clear varname fileout ind bobcall urlbase xmin xmax ymin ymax tpos tmin tmax xmin xmax ymin ymax

% SST and bathy
% save -mat get_stewart_SSTbathy_83048.mat %slat slon sst sst_3D stime %nt nx ny

% %% Wind
% 
% % use xtracto_3D_bdap. Zonal and meridional look like they call the same datasetname TQNux101day: they are already blended
% [windZ wlonZ wlatZ wtimeZ] = xtracto_3D_bdap([xmin xmax], [ymin ymax], [tmin tmax], 'qnux101day'); %25km zonal
% % [windM wlonM wlatM wtimeM] = xtracto_3D_bdap([xmin xmax], [ymin ymax], [tmin tmax], 'qnuy101day'); %25km meridional. 
% 
% wind_3D = squeeze(windZ); 
% [nt ny nx] = size(wind_3D);
% 
% save get_stewartOutput_wind_83051.mat nt nx ny windZ wlonZ wlatZ wtimeZ wind_3D  




