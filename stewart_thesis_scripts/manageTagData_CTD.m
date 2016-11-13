% manageTagData_CTD.m
%  this is called from compareTags.m
%
% This will import CTD data. User specifies which combo of data from:
%         MBARI ROV
%         NOAA Juvenile Rockfish

% This will also calculate horizontal distances, at the bottom of this
% script
%
% Outside Functions Called:
% plotCTD_Tags_2009 % J. Stewart
% plotCTD_Tags_Area % J. Stewart
% importCalCOFIbottles % J. Stewart
% horizontalDist_Tags % J. Stewart: calculate the minimum distance traveled
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 30-Mar-2010 13:10:22
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : manageTagData_CTD.m

reload = 1;
TagTempOverlay = 0;

Mfilename = mfilename;
dirCTD = '/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/_CTDs/';

%% Load CTD data: only run this cell if you want the data loaded.
% Data from ROVCTD and from JuvRockfish CTDs.

if reload
    disp(' ');
    castsROVmidwater = input('Import midwater MBARI ROVCTD data? (y=1, n=0) ');
    castsROVrect = 0;%input('Import rectangle MBARI ROVCTD data? (y=1, n=0) ');
    castsROVfall = input('Import tiburon fall MBARI ROVCTD data? (y=1, n=0) ');
    castsJuv = 0;%input('Import NOAA Juv Rockfish CTD data? (y=1, n=0) ');
    castsNov2010 = 0;%input('Import Gilly Lab Nov 2010 CTD data? (y=1, n=0) ');
    castsCalCOFI = input('Import CalCOFIbottle CTD data? (y=1, n=0) ');
    
    if ~castsCalCOFI
        plotAreaOrYear = input('Plot casts by area (1) or in 2009 (0)? ');
        disp(' ');
        disp('1: Monterey Tagging Area')
        disp('2: Monterey-Mexico tag')
        disp('3: Whole Juvenile Rockfish Coverage')
        disp('4: Big Sur')
        disp('5: Pismo Beach')
        disp('6: CalCOFI [Pismo Beach to Mx]')
        disp('7: SoCal Bight')
        disp('8: San Diego')
        Spatial = input('Which area would you like plotted? ');    % someday, make map?
%         colorData = input('Would you like to color by year (1) or by distance from shore (2)? ');
        
        % set up matrices, will concatenate them at the end, whether they are full or not
        matrixROVCTDup = [];
        matrixROVCTDr = [];
        matrixROVCTDf = [];
        matrixNov2010CTDdown = [];
        matrixJuvCTD = [];
        castIDROV = [];
        castIDROVrect = [];
        castIDROVfall = [];
        castIDJuv = [];
        castIDNov2010 = [];
        
        datasetIDtrix = [];
        if castsROVmidwater % load ROVCTD data
            load -mat Dg_ROVCTD_setup.mat % made in b_setupROVCTD.m by J. Stewart
            matrixROVCTD = matrixROVCTD(:,1:8);
            divenum = matrixROVCTD(:,1);
            unidivenum = unique(divenum);
            % clean up some ugly casts (ID'd in plotCTD_Tags_Area.m)
            badlog = divenum == 91.2 | divenum == 1791.1 | divenum == 1832.1 | divenum == 1833.1 | divenum == 1845.1;
            matrixROVCTD(badlog,:) = [];
            divenum = matrixROVCTD(:,1); % must do it again, because now missing dives.
            unidivenum = unique(divenum);
            
            for i = 1:length(unidivenum) % change ROV casts to just be upcast data. % TALK TO ROB/BRUCE ABOUT THIS
                divenuml = divenum == unidivenum(i);
                ToPlot = matrixROVCTD(divenuml,:);
                dep = ToPlot(:,5); % depth
                maxdep = max(dep);
                maxdepind = find(maxdep == dep);
                matrixROVCTDup = [matrixROVCTDup; ToPlot(maxdepind(1):end,:)];
            end
            castIDROV = 'ROV';
        end
        
        if castsROVrect % load ROVCTD data
            load -mat Dg_ROVCTD_setup_rect.mat % made in importROVCTD_rect.m by J. Stewart
            matrixROVCTD = matrixROVCTD(:,1:8);
            divenum = matrixROVCTD(:,1);
            unidivenum = unique(divenum);
%             % clean up some ugly casts (ID'd in plotCTD_Tags_Area.m)
%             badlog = divenum == 91.2 | divenum == 1791.1 | divenum == 1832.1 | divenum == 1833.1 | divenum == 1845.1;
%             matrixROVCTD(badlog,:) = [];
            divenum = matrixROVCTD(:,1); % must do it again, because now missing dives.
            unidivenum = unique(divenum);
            
            for i = 1:length(unidivenum) % change ROV casts to just be upcast data. % TALK TO ROB/BRUCE ABOUT THIS
                divenuml = divenum == unidivenum(i);
                ToPlot = matrixROVCTD(divenuml,:);
                dep = ToPlot(:,5); % depth
                maxdep = max(dep);
                maxdepind = find(maxdep == dep);
                matrixROVCTDr = [matrixROVCTDup; ToPlot(maxdepind(1):end,:)];
            end
            castIDROVrect = 'ROVrect';
            matrixROVCTD = [];
        end
        
        if castsROVfall % load tiburon fall data
            load -mat Dg_ROVCTD_setup_fall.mat % made in importROVCTD_rect.m by J. Stewart
            matrixROVCTD = matrixROVCTD(:,1:8);
            divenum = matrixROVCTD(:,1);
            unidivenum = unique(divenum);
            % clean up some ugly casts (ID'd in plotCTD_Tags_Area.m)
            badlog = divenum == 91.2 | divenum == 1791.1 | divenum == 1832.1 | divenum == 1833.1 | divenum == 1845.1;
            matrixROVCTD(badlog,:) = [];
            divenum = matrixROVCTD(:,1); % must do it again, because now missing dives.
            unidivenum = unique(divenum);
            
            % troubleshoot--keep this because removed bad files and putthem in folders labeled uniaa and unibb.
            a = find(oxyg <0);
            aa = divenum(a);
            uniaa = unique(aa);
            b = find(oxyg < 3 & temper > 13);
            bb = divenum(b);
            unibb = unique(bb);
            c = find(oxyg < 2 & temper > 3 & temper < 5);
            cc = divenum(c); %$ fix this weirdo too see screeenshot
            unicc = unique(cc);
            figure; hold on
%             plot(oxyg, temper, '.')
%             for i = 1:length(uniaa)
%                 ilog = uniaa(i) == divenum;
%                 plot(oxyg(ilog), temper(ilog), '.w')
%             end
%             for i = 1:length(unibb)
%                 ilog = unibb(i) == divenum;
%                 plot(oxyg(ilog), temper(ilog), '.r')
%             end

for i = 1:length(unidivenum) % change ROV casts to just be upcast data. % TALK TO ROB/BRUCE ABOUT THIS
                divenuml = divenum == unidivenum(i);
                ToPlot = matrixROVCTD(divenuml,:);
                dep = ToPlot(:,5); % depth
                maxdep = max(dep);
                maxdepind = find(maxdep == dep);
                matrixROVCTDt = [matrixROVCTDup; ToPlot(maxdepind(1):end,:)];
            end
            castIDROVfall = 'ROVfalltib';
            matrixROVCTD = [];
        end
        
        if castsJuv % load JuvRockfish data
            cd '/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/_CTDs/_JuvenileRockfishCruises'
            load -mat JuvCTD_2009.mat% JuvCTD_2008_2009.mat % made in importJuvCTD.m by J. Stewart
            castIDJuv = 'Juv';
        end
        
        if castsNov2010 % load Nov2010 data
            cd '/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/_CTDs/_FulmarVIPtrip_2010_11_30'
            load -mat Nov2010CTD.mat % made in import2010CTD.m by J. Stewart
            matrixNov2010CTD = matrixNov2010CTD(:,1:8);
            divenum = matrixNov2010CTD(:,1);
            unidivenum = unique(divenum);
            matrixNov2010CTDdown = [];
            for i = 1:length(unidivenum) % change casts to just be upcast data
                divenuml = divenum == unidivenum(i);
                ToPlot = matrixNov2010CTD(divenuml,:);
                dep = ToPlot(:,5); % depth
                maxdep = max(dep);
                maxdepind = find(maxdep == dep);
                matrixNov2010CTDdown = [matrixNov2010CTDdown; ToPlot(1:maxdepind(1),:)]; % do downcasts: the upcasts are really weird
            end
            castIDNov2010 = 'Nov2010';
        end
        
        % combine imported data matrices into one supermatrix
        matrixCTD = [matrixROVCTDup; matrixROVCTDr; matrixROVCTDf; matrixNov2010CTDdown; matrixJuvCTD];
        castID = [castIDROV castIDROVrect castIDROVfall castIDJuv castIDNov2010];
        
        % clean up  oxygen outliers
        divenum = matrixCTD(:,1);
        oxygCTD = matrixCTD(:,7);
        unidivenum = unique(divenum);
        
        oxygneglog = oxygCTD <= 0;
        oxygnegOxyg = oxygCTD(oxygneglog);
        oxygnegDiveNum = divenum(oxygneglog);
        oxygnegDiveNumUni = unique(oxygnegDiveNum); % perhaps too many to just throw out?
        
        for i = 1:length(oxygnegDiveNumUni) % remove any casts with oxygen problems
            ilog = oxygnegDiveNumUni(i) == divenum;
            matrixCTD(ilog,:) = NaN;
        end
        sum(matrixCTD(:,7) < 0);
        
        % simplify dataset
        divenum = matrixCTD(:,1);
        dateCTD = matrixCTD(:,2);
        latCTD = matrixCTD(:,3);
        lonCTD = matrixCTD(:,4);
        temperCTD = matrixCTD(:,6);
        oxygCTD = matrixCTD(:,7);
        unidivenum = unique(divenum(~isnan(divenum)));
%         clear oxyg temper lon lat dep depth divenuml dateCTDV matrixROVCTD
        
        % get info on individual dives
        dateCTDstart = [];
        for k = 1:length(unidivenum)
            templog = divenum == unidivenum(k);
            tempdate = dateCTD(templog);
            startdate = tempdate(1);
            dateCTDstart = [dateCTDstart; startdate];
        end
        dateCTDstartV = datevec(dateCTDstart);
        
    else % for CalCOFI
        plotAreaOrYear = 2;
        disp(' ');
        disp('1: Monterey Tagging Area')
        disp('2: Monterey-Mexico tag')
        disp('3: Whole Juvenile Rockfish Coverage')
        disp('4: Big Sur')
        disp('5: Pismo Beach')
        disp('6: CalCOFI [Pismo Beach to Mx]')
        disp('7: SoCal Bight')
        disp('8: San Diego')
        Spatial = input('Which area would you like plotted? ');    % someday, make map?
    end
end

%% Determine which CTDs to use, based on area/location

LocTagDeployCA = [ % CA tag locations
    %      -123.4831667, 37.91416667, ; % tagged:  83046_08
    -122.0263333, 36.70666667; % 83046_09
    -122.0556167, 36.58855; % 83052_09
    -122.0547667, 36.58445; % 64004_09
    -122.0569333, 36.5822; % 83051_09
    ];
LocTagPopUpCA = [ % CA tag locations
    %     -123.993, 37.645; % popoff:  83046_08
    -121.4955, 35.852; % 83046_09
    -123.2983333, 36.22888889; % 83052_09
    -122.0897222, 36.18277778; % 64004_09
    -118.4463889, 32.09194444 % 83051_09
    ];

xLimWest = -124.5;
xLimEast = -117; % unless there's a reason not to include land, can't think of one

switch Spatial  % identify area based on userinput in manageTagData_CTD % J. Stewart
    case {1} % Monterey Tagging Area
        yLimNorth = 37;
        yLimSouth = 35.5;
    case {2} % Monterey-Mexico tag % based on CalCOFI grid
        yLimNorth = 37;
        yLimSouth = 29.5;
    case {3} % Whole Juvenile Rockfish Coverage
        yLimNorth = 38.5;
        yLimSouth = 29.5;
    case {4} % Big Sur
        yLimNorth = 36;
        yLimSouth = 35.2;
    case {5} % Pismo Beach
        yLimNorth = 35.3;
        yLimSouth = 34.5;
    case {6} % CalCOFI [Pismo Beach to Mx]
        yLimNorth = 35.2;
        yLimSouth = 30;
    case {7} % SoCal Bight
        yLimNorth = 34.6;
        yLimSouth = 33.3;
    case {8} % San Diego
        yLimNorth = 33.3;
        yLimSouth = 30;
    otherwise
        error('Please choose a number corresponding to a location. ');
end
% to make boxes
SpatialLon = [xLimWest; xLimWest; xLimEast; xLimEast; xLimWest];
SpatialLat = [yLimSouth; yLimNorth; yLimNorth; yLimSouth; yLimSouth];
Locations = [SpatialLat SpatialLon];
laty = [min(SpatialLat) max(SpatialLat)];
lonx = [min(SpatialLon) max(SpatialLon)];

%% Plot CTDs based on Time (months around the 2009 tags) or Area

% setup for figures
maxd = 2000;
maxo = 8;
maxt = 20;
incrd = maxd*0.03;
incrt = maxt*0.03;

if plotAreaOrYear == 1
    AreaID = 'Area';
    plotCTD_Tags_Area
elseif plotAreaOrYear == 0
    AreaID = '2009tags';
    plotCTDuserinput % by J. Stewart: used for visualizing and choosing specific CTDs in region. Will also be called with OMLhistfig from manageTagData.m
    %     plotCTD_Tags_2009 % not sure whether this is necess
elseif plotAreaOrYear == 2
    AreaID = 'calCOFI';
    importCalCOFIbottles
    save -mat CalCOFIcasts.mat CalCOFIcasts LineNum StationNum UniLineNum UniStationNum
    calcofi_tags_ols  % J. Stewart: least squares analysis with tag 5. Extension of % calcofi_OLS.m, CalCOFI_Tag_OLS_dgf.m (D. Foley's version of CalCOFI_Tag_OLS.m)

    %     contourCalCOFIbottles % contours along the lines
%     plotCalCOFI_lines %Plot each line separately or temp, oxy, then temp v. oxy. and map. 
% plotCTD_TagOverlay % J. Stewart: overlay calcofi with tag. work from here with OLS

end

% horizontalDist_Tags % J. Stewart

%%

disp('Completed manageTagData_CTD.m')
% ===== EOF [manageTagData_CTD.m] ======
