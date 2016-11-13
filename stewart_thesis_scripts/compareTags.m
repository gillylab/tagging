% compareTags.m
%
% This loads .mat files made by importPDT_csv_JS.m and compares min/max
% temp/depth
%
% Outputs:
%
% Outside Functions Called:
% plotROVCTD_Tags.m
% makeLabels.m % A. Booth
% calcOMZtimes.m %J. Stewart
% calcDepthtimes.m %J. Stewart
% compareTagsFigs.m % J. Stewart
% compareTagsVertVel.m % J. Stewart
% compareTagsFigsBothLocs.m % J. Stewart
% compareTagsMBARI % J. Stewart
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 18-Mar-2010 14:10:42, 09-Oct-10 Revised to include GOC 09-Nov-10 Revised to include Temp,
% all CA and MX tags rerun on May 4 2011. Cleaned up 01-Aug-11 for all
% histograms
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : compareTags.m

Mfilename = mfilename;

%%
% plotROVCTD_Tags.m

TagsInCA = input('Comparison for California or Mexico? (1 = California, 0 = Mexico) ');
TagsVert = 1; % input('Load data with vertical migration (1) or without (0)? ');
DeploymentLength = input('run for entire deployment [calcDepthtimes] (1) or for first 3 days [calcOMZtimes] (0)? ');
figs = input('plot all compareTags figures? Y(1) or N(0): ');
dbarincr = input('depth bins by 25 or 50? (type value):'); % to compare with each other, do 25. to compare with MBARI, do 50
DL = 3; % what the hell is this? I changed it from where it was in calcDepthtimes because it made no sense. 

if TagsInCA
    dirCompare = '/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/CompareTagDir/_CA_compare/';
    % dirCompare = '/Users/jstewart/Dropbox/Thesis/Dg_Tagging/CompareTagDir/_CA_compare'
    rows = 2;
    cols = 3;
    tempMax = 20;
    locID = 'CA';
else
    dirCompare = '/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/CompareTagDir/_GOC_compare/';
    rows = 2;
    cols = 4;
    tempMax = 30;
    locID = 'GOC';
end

cd(dirCompare);
matfiles = dir(fullfile(dirCompare,'*.mat')); %find .mat files

%% Load .mat files made by importPDT_csv_JS.m

% import MinMaxTempDepth.mat files
if TagsVert
    matfilefun = @(x) char(regexp(x, '.+MinMaxTD.mat','match')); %files of interest
    vertID = '';
else
    matfilefun = @(x) char(regexp(x, '.+MinMaxTDnovert.mat','match')); %files of interest
    vertID = '-novert';
end

% read in histogram info
matfilename = cellfun(matfilefun,{matfiles(:).name}, 'UniformOutput',0); %find files of interest within cell array
matfiledex = cellfun(@isempty,matfilename);
matfilename = char(matfilename(~matfiledex));
matfilename = cellstr(matfilename);

CompareDT = [];
BoxPlotDepthtrix = [];
DepthHistoTrix = [];
BoxPlotTemptrix = [];
TempHistoTrix = [];
DepthHistotrixLength = [];
TempHistotrixLength = [];
DayCutoff = nan(length(matfilename), 1);
TagNameN = [];
for i = 1:length(matfilename)
    filetemp = ['load ' matfilename{i}];
    eval(filetemp);
    ii = matfilename{i};
    iii = str2double(ii(1));
    iv = str2double(ii(9:12));
    CompareDT = [CompareDT; Temptrix Depthtrix(:,3:4) repmat(iii,size(Temptrix),1)]; % Depthtrix is [DeployDates tagNumID Depthmin Depthmax] for each day of deployment
    TagNameN = [TagNameN; iii];
        
    % FOR BOXPLOTS
    BPtemp = Boxplotrix;
    if ~DeploymentLength % if run for first 3 days
        aa = BPtemp(:,7) == 1;
        bb = BPtemp(:,7) == 2;
        if sum(aa) < sum(bb)*0.3 % 30% seems like an appropriate length. less than 30% of a day, don't count it.
            DayCutoff(i) = 4;
        else
            DayCutoff(i) = 3;
        end
        alog = BPtemp(:,7) > DayCutoff(i); % deployment day index
        BPtemp(alog,:) = [];
    else
        if TagsInCA == 1 && i == 3 % CA tag 4: it floats on surface
            alog = BPtemp(:,7) > 3; % floats on day 4
            BPtemp(alog,:) = [];
        end
    end
    BoxPlotDepthtrix = [BoxPlotDepthtrix; BPtemp(:,1:4) repmat(iii,length(BPtemp),1) BPtemp(:,7)]; % these end bits are for identification!
    BoxPlotTemptrix = [BoxPlotTemptrix; [BPtemp(:,1:2) BPtemp(:,5:6)] repmat(iii, length(BPtemp), 1) BPtemp(:,7)]; % these end bits are for identification!
end  

% FOR HISTOGRAMS % minimal for histograms: need to redo individually.
% this works since BoxPlotDepthtrix is already tailored above for 3-day or not
dbar = 0:dbarincr:max([BoxPlotDepthtrix(:,3);BoxPlotDepthtrix(:,4)]); % 3 is night, 4 is day
tbar = (floor(min([BoxPlotTemptrix(:,3);BoxPlotTemptrix(:,4)])):0.5:max([BoxPlotTemptrix(:,3);BoxPlotTemptrix(:,4)])); % 5 is night, 6 is day

BoxPlottrixHeaders = ['TagNum', 'DEPLOYINDEXN', 'DepthN_OR_TempN', 'DepthD_OR_TempD', 'TagID', 'DayOfDeployment'];
HistoTrixHeaders = ['TagNum', 'DeployYear', 'tbarORdbar', 'percTN_OR_percDN', 'percTD_OR_percTN'];
TagIDs = unique(BoxPlotTemptrix(:,5));

% read in raw data for other analyses
matfilefun = @(x) char(regexp(x, '.+rawTD.mat','match')); %files of interest
matfilename = cellfun(matfilefun,{matfiles(:).name}, 'UniformOutput',0); %find files of interest within cell array
matfiledex = cellfun(@isempty,matfilename);
matfilename = char(matfilename(~matfiledex));
matfilename = cellstr(matfilename);

%FOR RAW DATA
Rawtrix = [];
for i = 1:length(matfilename)
    Rtrix = []; % this happens here because it didn't in processSeriesData.m
    filetemp = ['load ' matfilename{i}];
    eval(filetemp);
    ii = matfilename{i};
    iii = str2double(ii(1));
    
    Rtrix = [Rtrix; DateRaw DepthRaw TempRaw TagRaw DeployDayRaw repmat(iii,length(DateRaw),1)]; % to label
    
    if ~DeploymentLength % if run for first 3 days
        aa = Rtrix(:,5) == 1;
        bb = Rtrix(:,5) == 2;
        if sum(aa) < sum(bb)*0.3 % 30% seems like an appropriate length. less than 30% of a day, don't count it.
            DayCutoff(i) = 4;
        else
            DayCutoff(i) = 3;
        end
        alog = Rtrix(:,5) > DayCutoff(i); % deployment day index
        Rtrix(alog,:) = [];
    else
        if TagsInCA == 1 && i == 3 % CA tag 4: it floats on surface. maybe even cut it earlier, because dead?
            alog = Rtrix(:,5) > 3; % floats on day 4
            Rtrix(alog,:) = [];
        end
    end
    Rawtrix = [Rawtrix; Rtrix];  
    
%     % troubleshoot jan 19 2012
%     figure
%     plot(TempRaw,DepthRaw, '.')
%     title(ii)
%     set(gca,'YDir', 'reverse')
end

if DeploymentLength % full deployment
    dirBoth = '/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/CompareTagDir/_CA_GOC/';
    labDeploy = [];
else % 3day
    dirBoth = '/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/CompareTagDir/_CA_GOC_3day/';
    labDeploy = '3day';
    dirCompare = ([dirCompare '/_' labDeploy]);
end

% %save as .mat so can make figures with both in compareTagsFigs.m. ONLY DO WITH FULL DEPLOYMENT, NOT 3day
% filenameSave = ['AllTagsCompareData' locID '.mat'];
% eval (['save -mat ' filenameSave ' TagIDs', ' Rawtrix', ' CompareDT', ' BoxPlotDepthtrix', ...
%      ' DepthHistoTrix', ' BoxPlotTemptrix', ' TempHistoTrix', ' DepthHistotrixLength', ' TempHistotrixLength']);

%% now:

L = 80;
if DeploymentLength
    calcDepthtimes
else
    calcOMZtimes % calculate times in OMZ to create something like Figure 12 in Gilly etal 2006
end

if figs
    compareTagsFigs% J. Stewart
    compareTagsVertVel % J. Stewart
    compareTagsFigsBothLocs% J. Stewart % THIS CAN RUN ON ITS OWN--Can clear all first
end

%%

disp('Completed compareTags.m')
% ===== EOF [compareTags.m] ====



