% habitatID_eachTag.m
%
% This is called by manageTagData.m
%
% this will load SST-bathy data for the tag specified by the user. Gives
% user the choice to prepare or run the model.
% either prepare for Fortran or run things post Fortran. See
% FortranHowTo.txt for more details.
%
% Outside Functions Called:
% determineShallowCutoff.m (J. Stewart)
% compare_satellite.m (J. Stewart)
% probHabitat.m (J. Stewart)
% mapNEPac (A. Booth)
% radiusApproach_Tags (J. Stewart)
% randomWalk_Tags (J. Stewart)
% mappingPostFortran (J. Stewart)
% mappingPostFortran_dist (J. Stewart)
%
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 11-Aug-2010 12:01:35, edited 07-Mar-2011 after working with D. Foley: can't use polyval, use quantiles. Older code saved in habitatID.m.
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : habitatID_eachTag.m

dirHinke = '/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/_Hinke/';
dirFortran = '/Users/juliastewart/fortran/';

%% user input to identify tags

disp(' ');
disp('You will select a tag in order to run the Hinke analysis.');
disp('64004: (Nov Monterey Tag)');
disp('83046: (Oct-Nov Monterey Tag)');
disp('83048: (Dec Monterey/Maui Tag)');
disp('830462008: (2008 Cordell Bank Tag) [no space]');
disp('83051: (Nov Monterey-Mexico Tag)');
disp('83052: (Nov Monterey Tag)');
TagHinke = input('===> Please identify the tag number you would like to work with from above: ');

% input about what to run, whether to make figures
prep = input('==> prepare for fortran (1) or run model post Fortran? (0): ');
if ~prep
    mm = input('===> make probability maps, post Fortran? 1(Y) or 0(N): ');
    mmm = input('===> make distance figures, post Fortran? 1(Y) or 0(N): ');
    mmmm = input('===> create file for ArcGIS, post Fortran? 1(Y) or 0(N): ');
end

%% Work with user input

% DailyShallow files created in tagDailyHistos_Shallow.m and tagDailyHistos_ShallowPDT.m DON"T USE determineShallowCutoff ANYMORE!!!
% JUST READ IN Shallowest in below. It works with the
% get_stewart.m created these get_stewartOutput files: uses xtractomatic to access SST data

n = 5;
if TagHinke == 83046 % do this so will read in correct 83046 file
    cd('/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/_CA_Tags/83046_09');
    filenameLoad = ['DailyShallowest' num2str(n) '_' num2str(TagHinke) '.mat']; % load data depending on the cutoff you want
elseif TagHinke == 830462008
    cd('/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/_CA_Tags/83046_08');
    filenameLoad = ['DailyShallowest' num2str(n) '_' num2str(83046) '.mat']; % load data depending on the cutoff you want
else
    filenameLoad = ['DailyShallowest' num2str(n) '_' num2str(TagHinke) '.mat']; % load data depending on the cutoff you want
end

eval (['load -mat ' filenameLoad])
cd(dirHinke);

switch TagHinke
    
    case 64004 % tag 4
%         load -mat get_stewart_SST_wholearea.mat % just made this temporarily to make a sst map with patrick
        load -mat get_stewart_SSTbathy_64004.mat % was get_stewartOutput_64004.mat old way without bathy
        LocTagCA_Hinke = [-122.0547667, 36.58445; -122.0897222, 36.18277778];
        Deployment_day = 3; % start date Nov 6
        array3popup = [Deployment_day 1]; % this works with Deployment_day since there is actually 1 hour on day 4, so can give 1 last location then.
        HinkeDates = datenum('6-Nov-2009'):1:datenum('6-Nov-2009')+Deployment_day;% start date Nov 6
        NEPacBuffer = 1.5;
        dirf = '/Users/juliastewart/fortran/64004/';
        
    case 83046 % tag 2
        load -mat get_stewart_SSTbathy_83046.mat
        LocTagCA_Hinke = [-122.0263333, 36.70666667; -121.4955, 35.852];
        Deployment_day = 12.7; % start date Sept 30
        array3popup = [floor(Deployment_day) 1]; % these should all be -1 than deployment day so that it works with Day 0 as the deployment date
        HinkeDates = datenum('30-Sep-2009'):1:datenum('30-Sep-2009')+Deployment_day;% start date Nov 6
        NEPacBuffer = 1.5;
        dirf = '/Users/juliastewart/fortran/83046_09/';
        
    case 83048 % tag 2
        load -mat get_stewart_SSTbathy_83048.mat
        LocTagCA_Hinke = [-122.0455667, 36.69466667; -122.0455667, 36.69466667]; % no popup data here, put same for both
        Deployment_day = 6; % start date Dec 5
        array3popup = [floor(Deployment_day) 1]; % these should all be -1 than deployment day so that it works with Day 0 as the deployment date
        HinkeDates = datenum('05-Dec-2009'):1:datenum('11-Dec-2009')+Deployment_day;% start date Nov 6
        NEPacBuffer = 1.5;
        dirf = '/Users/juliastewart/fortran/83048_09/';
        
    case 83051 % tag 5
        load -mat get_stewart_SSTbathy_83051.mat
        %         load -mat calcofi_tags_olsAscent.mat % made in calcofi_tags_ols.m: least squares output for tag 83051
        load -mat calcofi_tags_olsA.mat
        load -mat calcofi_tags_olsD.mat
        %         load -mat calcofi_tags_olsDescent.mat % made in calcofi_tags_ols.m: least squares output for tag 83051
        ccA = [1 0.2 0;... % ascent colors
            0.87 0.49 0; ...
            0.87 0.49 0; ...
            1 0.6 0];
        ccD = [0.6 0 0.4; ... % descent colors
            1 0 1;
            1 0.2 1; ...
            1 0.4 1];
        LocTagCA_Hinke = [-122.0569333, 36.5822; -118.4463889, 32.09194444];
        Deployment_day = 17.6;
        array3popup = [floor(Deployment_day) 1]; % these should all be -1 than deployment day so that it works with Day 0 as the deployment date
        HinkeDates = datenum('6-Nov-2009'):1:datenum('6-Nov-2009')+Deployment_day;% start date Nov 6
        NEPacBuffer = 0.15;
        dirf = '/Users/juliastewart/fortran/83051/';
        
    case 83052 % tag 3
        load -mat get_stewart_SSTbathy_83052.mat
        LocTagCA_Hinke = [-122.0556167, 36.58855; -123.2983333, 36.22888889];
        Deployment_day = 5.7; % start date Nov 6
        array3popup = [floor(Deployment_day) 1]; % these should all be -1 than deployment day so that it works with Day 0 as the deployment date
        HinkeDates = [datenum('6-Nov-2009'):1:datenum('6-Nov-2009')+Deployment_day datenum('6-Nov-2009')+Deployment_day+1];% start date Nov 6
        NEPacBuffer = 1.5;
        dirf = '/Users/juliastewart/fortran/83052/';
        
    case 830462008 % tag 1
        load -mat get_stewart_SSTbathy_83046_08b.mat
        LocTagCA_Hinke = [-123.4831667, 37.91416667; -123.993, 37.645];
        Deployment_day = 2.7; % start date Oct 12
        array3popup = [floor(Deployment_day) 1]; % these should all be -1 than deployment day so that it works with Day 0 as the deployment date
        HinkeDates = datenum('12-Oct-2009'):1:datenum('12-Oct-2009')+Deployment_day;% start date Nov 6
        NEPacBuffer = 1.5;
        dirf = '/Users/juliastewart/fortran/83046_08/';
        
        
    otherwise
        error('Please identify a tag number from the list above');
end

divemax = DepthTagMax;

%% Location information for mapNEPac.m. There is a better way to index this. Do it later.

LocTagDeployCA = [ % CA tag locations
    -123.4831667, 37.91416667, ; % tagged:  83046_08
    -122.0263333, 36.70666667; % 83046_09
    -122.0556167, 36.58855; % 83052_09
    -122.0547667, 36.58445; % 64004_09
    %     -122.0569333, 36.5822; % 83051_09
    -122.0455667, 36.69466667; % 83048_09
    ];

LocTagPopUpCA = [ % CA tag locations
    -123.993, 37.645; % popoff:  83046_08
    -121.4955, 35.852; % 83046_09
    -123.2983333, 36.22888889; % 83052_09
    -122.0897222, 36.18277778; % 64004_09
    %     -118.4463889, 32.09194444 % 83051_09
    -122.0455667, 36.69466667; % 83048
    ];

if TagHinke == 83051
    LocTagDeployCA = [LocTagDeployCA; -122.0569333, 36.5822]; % 83051_09 deploy
    LocTagPopUpCA = [LocTagPopUpCA; -118.4463889, 32.09194444]; % 83051_09 popup
end

LocationsCA = [LocTagDeployCA; LocTagPopUpCA];
Locations = [LocationsCA(:,2) LocationsCA(:,1)];

%% make SST maps for each day of each tag, based on quantiles. From D. Foley, 07-Mar 2011

% just take the range, not quantile
qtrix = [];
dtrix = []; % deepest depth trix
for i = 1:size(DepthTempTagDay,2)/2 % this should be the same length as stime
    qmin = min(DepthTempTagDay(:,i*2));
    qmax = max(DepthTempTagDay(:,i*2));
    qtrix = [qtrix; qmin qmax];
    dmin = min(DepthTempTagDay(:,i*2-1));
    dmax = max(DepthTempTagDay(:,i*2-1));
    dtrix = [dtrix; dmin dmax];
end

if TagHinke == 83046
    qtrix(1,:) = [12.5 12.7]; % enter CTD data for that first super deep value
    dtrix(1,:) = [0 20];
end

disp('see CompareShallowCutoffs.xls for some comparisons of shallowness');
var(dtrix,0,2)
var(qtrix,0,2)

std(dtrix,0,2)
std(qtrix,0,2)


%% prepare or run the model

if prep % prepare
    
%     radiusApproach_RicardoPrep % J. Stewart 24 May 2011
    radiusApproach_RicardoPrepNoBathy % J. Stewart 10 Apr 2011
    
else % run the model
    if TagHinke == 83048
        radiusApproach_Tags % map it without probabilities
    else
        
        cd(dirf); % read output from fortran script
        
        % map data
        filename = 'SquidSST_meshb.txt';
        fid=fopen(filename);
        mapData = textscan(fid, '%n %n %n','HeaderLines', 1, 'Delimiter','\t');
        fclose(fid);
        mapLat = mapData{1,1};
        mapLon = mapData{1,2};
        mapOcean = mapData{1,3};
        msk = mapOcean==1;
        mapLat = mapLat(msk);
        mapLon = mapLon(msk);
        
        % probability data
        files = dir([dirf '/*_positions.txt']);
        filename = files.name;
        fid=fopen(filename);
        probData = textscan(fid, '%n %n %n %n','HeaderLines', 1, 'Delimiter','\t');
        fclose(fid);
        probDay = probData{1,1};
        probLat = probData{1,2};
        probLon = probData{1,3};
        probProb = probData{1,4};
        probDayUni = unique(probDay);
        
        
        swimmaxS = regexp(filename, '_', 'split');
        swimmax = str2double(swimmaxS(2));
        
        if mm % make probability maps
            mappingPostFortran % calls manageTagData_CTD if tag 83051 % radiusApproach_Tags % has been replaced with mappingPostFortran.m % J. Stewart
        end
        
        if mmm % make figures from probabilty calculations
            mappingPostFortran_dist % calcs dist between each day and distoffshore
        end
        
        if mmmm % make files for ArcGIS
            createTagFilesGIS
            disp(' ');
            disp(' ');
            disp('---- must delete spaces and change \t to , and save as .csv in TextWrangler ---- ');
            disp(' ');
            disp(' ');
        end
    end
end


%%

disp('Completed habitatID_eachTag.m')
% ===== EOF [habitatID_eachTag.m] ======
