% importTagData_csv.m
%
% This program will read in tag data and will call 2 programs to process it, depending on whether it is
% timeseries data (ProcessSeriesData.m) or PDT data (ProcessPDTdata.m) Currently 83046_09 is the only PDT tag. 
% The user will be asked to specify which type of data the tag has. 
%
% importTagData_csv.m (this program is based on importSeries.m by A. Booth but now to
% imports .csv files exported by Wildlife Computers software>Argos Processing>"Argos Message Decoder" (from WC-DAP) instead of .xls files
% 
% Outside functions called: 
%   subsampleArchivalTag.m % J. Stewart
%   ProcessPDTdata % J. Stewart and A. Booth
%   ProcessSeriesData % J. Stewart and A. Booth
%
% Created Oct 17th 2008 from importSeries.m; A. Booth, ashley@boothfamily.com
% Changed to importTagData_csv 27-Oct-2009 J. Stewart jules32@gmail.com.  
% This version processes timeseries, whether recovered or downloaded, or downloaded PDT data. 
% The previous/old version saved as importTagData_csvold.m

Mfilename = mfilename;

disp(' ');
disp('===> You will be processing data from a time series tag.');
TagLocation = input('Where was this tag deployed? (1 = California, 0 = GOC, 2 = Magdalena Bay) ');
TagRecovered = input('Was this tag recovered? Y=1 (1-second data), N=0 (75-second data):  ');
OMLhistfig = 0;%input('===> Would you like to plot depth histos with OML? (Y=1, N=0) '); 
ShallowCutoff = 20; % used in SST, Hinke calculations % 50 would be as deep as you could go, depending on the winds. See stats below on each cutoff. You lose quite a bit of the data the shallower you are, and temps dont' change much. 

if TagLocation == 1
   % cd('/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/_CA_Tags/');
    cd('/Users/jstewart/Dropbox/Thesis/Dg_Tagging/_CA_Tags/');
    CompareTagDir = '/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/CompareTagDir/_CA_compare/';
    colorsLoc = [4 16]; % California temps 
elseif TagLocation == 0
    % cd('/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/_GOC_Tags/');
    cd('/Users/jstewart/Dropbox/Thesis/Dg_Tagging/_GOC_Tags/');
    CompareTagDir = '/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/CompareTagDir/_GOC_compare/';
    disp('-->for tags with deployments crossing months (labeled blue), you must modify processSeriesData.m and plotSeriesData.m slightly');
    disp('search for {== DeployDays(i)} and replace with {== i}');
    userOK = input('press any key to acknowledge this');
    colorsLoc = [5 35]; % Mexico temps
elseif TagLocation == 2
    cd('/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/_Magdalena_Tags/');
    CompareTagDir = '/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/CompareTagDir/_Magdalena_compare/';
end
    
%% Import time series data, either from recovered or uploaded tag

if TagRecovered
    oldTag = input('was this tag deployed before 2009? 1(Y), 0(N): ');
    if oldTag
        [Depth_unsmooth,TempTag,Light,Seriesdate,time,dir1,file,name] = importArchFile;% (A. Booth 2009) % this will make you crop at popoff
        else
        [Depth_unsmooth,TempTag,Light,Seriesdate,time,dir1,file,name] = importArchFile_v2;% (modified J. Stewart 2011 for new tag formats)
    end
    if TagLocation == 1
        dir1 = [dir1 '/'];
    end
    % Smooth the data and get vertical velocities
    Depth = sgolayfilt(Depth_unsmooth,6,33); %smooths quantized trace while keeping large peaks
    DepthTagtemp = sgolayfilt(Depth,3,43);%sgolayfilt(Depth,2,13); %apply filter again to smooth little wiggles
    DepthTag = DepthTagtemp*-1;
    Vertvel = diff(Depth)./diff(time);%Vertvel = diff(Depth,1); %1st derivative
    tagdex = regexp(name, '(\d{5})');
    tagNumS = name(tagdex:tagdex+4); % this just works for 5-numbered PTT tags
    
    % import -Argos to find popoff information to use for sunrise and sunset
    csvfiles = dir(fullfile(dir1,'*.csv')); %find .csv files    
    Argofilefun = @(x) char(regexp(x, '.+-Argos.csv','match')); %files of interest
    Argofilename = cellfun(Argofilefun,{csvfiles(:).name}, 'UniformOutput',0); %find files of interest within cell array
    Argofiledex = cellfun(@isempty,Argofilename);
    Argofilename = char(Argofilename(~Argofiledex));
    
    % account for GOC tags whose software was pre-.csv files
    if length(Argofilename) > 1 % if there was such an Argofilename.csv
        fid=fopen(fullfile(dir1,Argofilename));
        line = 1;
        while line<5 %get headers
            tline = fgetl(fid);
            if ~isempty(tline) && strcmp(tline(1),'D'), break, end %get header names
            line = line+1;
        end
        fclose(fid);
        Tline = ['''' regexprep(tline,',',''';''') '''']; %reformat headers
        Argos.headernames = eval(['{' Tline '}']); %make into string array
        Argos.headernames = Argos.headernames';
        numcol = length(Argos.headernames);
        
        %import column data
        fmat = '%s  %n  %s	%s	%n	%n	%n	%n	%n	%n	%s  %s'; %inicial format of each columm
        for s = 1:numcol-12 %make colheaderFormat line adding on %n as long as headers are
            fmat = strcat(fmat,' %n');
        end
        fid=fopen(fullfile(dir1,Argofilename));
        Argos.data = textscan(fid, fmat,'HeaderLines', line, 'Delimiter',',');
        fclose(fid);
        
        datecol = find(cellfun(@(x) strcmp(x,'Date'),Argos.headernames)==1); %find Date column
        Argos.data{datecol} = double(Argos.data{datecol});
        latcol = find(cellfun(@(x) strcmp(x,'Latitude'),Argos.headernames)==1); %find Lat column
        latdata = Argos.data{1,latcol};
        lat = latdata(1); %getting the first lat (=popup)
        loncol = find(cellfun(@(x) strcmp(x,'Longitude'),Argos.headernames)==1); %find lon column
        londata = Argos.data{1,loncol};
        lon = londata(1); %getting the first lon (=popup)
    else
        matfiles = dir(fullfile(dir1,'*.mat')); %find .mat files
        PopUpfilefun = @(x) char(regexp(x, '.+tagloc.mat','match')); %files of interest: location mat file
        PopUpfilename = cellfun(PopUpfilefun,{matfiles(:).name}, 'UniformOutput',0); %find files of interest within cell array
        PopUpfiledex = cellfun(@isempty,PopUpfilename);
        PopUpfilename = char(PopUpfilename(~PopUpfiledex));
        filepresent = exist([dir1 PopUpfilename]);
        if filepresent == 2
            load([dir1 PopUpfilename], '-mat')
            lat = tagloc(1);
            lon = tagloc(2);
        else
            disp('---->> no popup location data for this tag, using CA deployment info')
            lat = 36.69; % these are deployment positions for Tag 6: 83048_09
            lon = -122.05;
        end
    end
    
else % tag not recovered
    TagTimeSeries = input('===> Is there time series data downloaded for this tag? Y=1, N=0 ');
    
    if TagTimeSeries % unrecovered tag with downloaded time series data
        disp(' ')
        disp('==> Select folder where *-Series.csv file is contained')
        dir1 = uigetdir;
        dir1 = [dir1 '/'];
        cd(dir1)
        
        % find files
        csvfiles = dir(fullfile(dir1,'*.csv')); %file .csv files
        Seriesfilefun = @(x) char(regexp(x, '.+-Series.csv','match')); %files of interest
        Seriesfilename = cellfun(Seriesfilefun,{csvfiles(:).name}, 'UniformOutput',0); %find files of interest within cell array
        Seriesfiledex = cellfun(@isempty,Seriesfilename);
        Seriesfilename = char(Seriesfilename(~Seriesfiledex));
        
        % import -Series.csv
        fid=fopen(fullfile(dir1,Seriesfilename));
        line = 1;
        while line<5 %get header names
            tline = fgetl(fid);
            if ~isempty(tline) && strcmp(tline(1),'D'), break, end
            line = line+1;
        end
        fclose(fid);
        Tline = ['''' regexprep(tline,',',''';''') '''']; %reformat
        Series.headernames = eval(['{' Tline '}']); %make into string array
        numcol = length(Series.headernames);
        
        %import column data
        fmat = '%s %n %s %s %s %n %n %n %n %n %n %n';
        for s = 1:numcol-7 %make colheaderFormat line adding on %n as long as headers are
            fmat = strcat(fmat,' %n');
        end
        fid=fopen(fullfile(dir1,Seriesfilename));
        Series.rawdata = textscan(fid, fmat,'HeaderLines', line, 'Delimiter',',');
        fclose(fid);
        
        %find pop line; get lat/lon of popup (works best if the Minimum Argos Location Class is better than 0 or 1-3)
        loccol = find(cellfun(@(x) strcmp(x,'LocationQuality'),Series.headernames)==1); %find LocationQuality column
        locqual = ~isnan(Series.rawdata{1,loccol}); % series has nans %locqual = cell2mat(cellfun(@(x) ~isempty(x),Series.rawdata{1,loccol},'UniformOutput',0));% %make sure locqual is numeric
        locqualsum = sum(locqual);
        if locqualsum > 0 % when it did pop-up
            locqual2 = find(locqual==1);
            Seriespopline = locqual2(1)-1; %end of deployment before start transmitting
            Series.rawdatacrop = cellfun(@(x) x(1:Seriespopline),Series.rawdata,'UniformOutput',0); %crop at pop off
            latcol = find(cellfun(@(x) strcmp(x,'Latitude'),Series.headernames)==1); %find Lat column
            latdata = Series.rawdata{1,latcol};
            lat = latdata(Seriespopline+1); %getting the first lat (=popup)
            loncol = find(cellfun(@(x) strcmp(x,'Longitude'),Series.headernames)==1); %find lon column
            londata = Series.rawdata{1,loncol};
            lon = londata(Seriespopline+1); %getting the first lon (=popup)
            
        else % when no pop-up data
            Series.rawdatacrop = Series.rawdata;
            lat = 36.69; % these are deployment positions for Tag 6: 83048_09
            lon = -122.05;
        end
        
        % date stuff
        daycol = find(cellfun(@(x) strcmp(x,'Day'),Series.headernames)==1); % find Day column
        timecol = find(cellfun(@(x) strcmp(x,'Time'),Series.headernames)==1); % find Time column
        dayvec = datevec(Series.rawdatacrop{1,daycol});
        timevec = datevec(Series.rawdatacrop{1,timecol});
        Seriesdate = datenum([dayvec(:,1:3), timevec(:,4:6)]);
         
        %get PTT number as a string
        pttcol = find(cellfun(@(x) strcmp(x,'Ptt'),Series.headernames)==1); %find ptt column
        pttdata = Series.rawdata{1,pttcol};
        tagNumS = num2str(pttdata(1)); %char(regexp(Seriesfilename, '\d{5}(?=-\w+.csv)','match'));
        
        depthcol = find(cellfun(@(x) strcmp(x,'Depth'),Series.headernames)==1); %find Depth column
        DepthTag = Series.rawdatacrop{1,depthcol};
        tempcol = find(cellfun(@(x) strcmp(x,'Temperature'),Series.headernames)==1); %find Depth column
        TempTag = Series.rawdatacrop{1,tempcol};
        
        Vertvel = diff(DepthTag)./diff(Seriesdate);%Vertvel = diff(Depth,1); %1st derivative
        
    else % unrecovered tag with downloaded PDT data (without downloaded time series data)
       processPDTData % by J. Stewart   
    end 
end

%% Process Time Series Data

pf = input('Make all figures? 1=y/0=n: ');
sampleint = 75; % timeseries sampling rate. identify here for later use in calcOMZtimes.m, and in subsampleArchivalTag.m

if TagRecovered || TagTimeSeries
    TagDataType = 1;%input('Would you like nights defined from dusk to dawn (1) or as surface (0)?');
    if TagRecovered 
       subsampleArchivalTag % J. Stewart % this assumes 1 second sampling frequency
    end
    processSeriesData % by J. Stewart
    saveTagDataAsMat % by J. Stewart, this will make a mat file that will run with Ashley's data
end

if ~exist('Light', 'var')
    Light = nan(length(DepthTag), 1);
end

% trackit = [SeriesdateV DepthTag Light TempTag];
% trackitName = [tagNameMat '_trackit.txt'];
% dlmwrite(trackitName, trackit, '\t');

%% save variables for tag comparison

tagNumD = str2num(tagNumS);
Temptrix = [DeployDays repmat(tagNumD, length(DeployDates), 1) Temptrix];
Depthtrix = [DeployDays repmat(tagNumD, length(DeployDates), 1) Depthtrix];
% 
% save([CompareTagDir tagName '_MinMaxTD.mat'] , 'Temptrix', 'Depthtrix', 'Boxplotrix', 'DepthHistotrix', 'TempHistotrix') % or '_MinMaxTDnovert.mat'

DateRaw = Seriesdate;
DepthRaw = DepthTag;
TempRaw = TempTag;
TagRaw = repmat(tagNumD,length(DateRaw),1);
DeployDayRaw = DEPLOYINDEX-DeployDays(1)+1;
% save([CompareTagDir tagName '_rawTD.mat'], 'DateRaw', 'DepthRaw','TempRaw', 'TagRaw', 'DEPLOYINDEXN', 'tagNumS', 'DeployDayRaw' )

% save([CompareTagDir tagName '_vertMigrations.mat'], 'DepthDescenti', 'TempDescenti', 'SeriesdateNtimesiD', 'DepthAscenti', 'TempAscenti', 'SeriesdateNtimesiA', 'tagNumS')
% save([CompareTagDir tagName '_vertMigrationsIndivid.mat'], 'AscentDatetrixVtimes', 'AscentDepthtrix', 'DescentDatetrixVtimes', 'DescentDepthtrix', 'tagNumS') This is only if you're running the Individ_June2011 version of process and plotSeriesData

%%

disp('Completed importTagData_csv.m') 

