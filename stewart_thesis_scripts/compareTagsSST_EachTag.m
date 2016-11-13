% compareTagsSST_EachTag.m
%  
% this is called by manageTagData.m  
%  
% Cell titles:
%  
% Outside Functions Called: 
%  
% AUTHOR    : J. Stewart jules32@gmail.com 
% DATE      : 22-Feb-2011 18:56:36  
% Revision  : 1.00  
% DEVELOPED : 7.9.0.529 (R2009b) OSX 
% FILENAME  : compareTagsSST_EachTag.m 

% compareTagsSST.m

Mfilename = mfilename;

TempBinIncr = .5; % probably smallest you'd want to do is .5, check T error on tags
DepthBinIncr = 25;

%%%
ShallowCutoff = 30; % 50 would be as dep as you could go, depending on the winds. See stats below on each cutoff. You lose quite a bit of the data the shallower you are, and temps dont' change much. 
%%% 

%% load data from all tags

dirCompare = '/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/CompareTagDir/_CA_compare'; 
cd(dirCompare);

% import rawTD.mat files 
matfiles = dir(fullfile(dirCompare,'*rawTD.mat')); %find .mat files % CURRENTLY DOES NOT INCLUDE 2008 TAG (moved to different folder)

% load in each file and save data in one big matrix. 
DepthAll = []; % all Depths
TempAll = []; % all Temps
DatesAll = []; % all Dates
DeployDatesAll = []; % deployment dates, but actually just each individual day of tagging
TagNames = []; % individual tag names
DataLength = []; % data length for each tag. 
TagID = [];
for i = 1:length(matfiles)
    infile = matfiles(i).name;
    load(infile)
    DepthAll = [DepthAll; DepthRaw]; 
    TempAll = [TempAll; TempRaw];
    DatesAll = [DatesAll; DateRaw];
    DeployDatesAll = [DeployDatesAll; DeployDates];
    TagNames = [TagNames; TagRaw]; 
    DataLength = [DataLength; length(DepthRaw)];
    TagID = [TagID; repmat(i, length(DepthRaw), 1)]; % add i's so can index the tag properly (ie year)
end

TagNamesUni = unique(TagNames(:,1)); % this wouldn't distinguish between 83046_2008 and 83046_2009
TagIDsUni = unique(TagID);
    
clear infile DeployDates DepthTag TempTag i 

%% make a filter for shallow depths

% for all tags
DepthShallowlog = DepthAll <= ShallowCutoff; % ShallowCutoff set above
DepthShallow = DepthAll(DepthShallowlog,1);
DatesAllShallow = DatesAll(DepthShallowlog,1);
TempShallow = TempAll(DepthShallowlog,1);

length(TempAll) % total amount of data counts for all tags
DataLengthShallow = length(TempShallow) % right now not for individual tags: but need to see that. ***TO DO***

% length(TempAll)     % with all 5 tags         29605
% ShallowCutoff = 30; length(TempShallow)       3909    13%
% ShallowCutoff = 40; length(TempShallow)       6758    23%
% ShallowCutoff = 20; length(TempShallow)       1807    06%
% ShallowCutoff = 10; length(TempShallow)       762     03%


% for each tag, make a matrix: the column number is the tagID.

% dtCell = {TagIDsUni, [1 2 3]} % not useful here now but
% maxVal = max(cellfun(@max,dtCell)); % this will tell you the maximum size of all of them. From stackoverflow.com 

% not the ideal way to do this, but for now...
for i = 1:length(TagIDsUni)
    ilog = TagID == TagIDsUni(i);
    iDepthShallow = ilog + DepthShallowlog;
    iDepthShallowlog = iDepthShallow == 2;
    DepthEach = DepthAll(iDepthShallowlog);
    DatesEach = DatesAll(iDepthShallowlog);
    TempEach = TempAll(iDepthShallowlog);
    if i == 1
        Tag1_Shallow = [DepthEach DatesEach TempEach]; % so the column order is Depth, Date, Temp. All at shallow cutoff. 
    elseif i == 2
        Tag2_Shallow = [DepthEach DatesEach TempEach];
    elseif i == 3
        Tag3_Shallow = [DepthEach DatesEach TempEach];
    elseif i == 4
        Tag4_Shallow = [DepthEach DatesEach TempEach];
    end
end
TagShallowNames = {'Tag1_Shallow', 'Tag2_Shallow', 'Tag3_Shallow', 'Tag4_Shallow'};

clear DepthEach DatesEach TempEach ilog iDepthShallow iDepthShallowlog


%% Histograms of Shallow. this is modified from processSeriesData.m

for i = 1%:length(TagShallowNames)
    DepthShallowA = eval(TagShallowNames{i});
    DepthShallow = DepthShallowA(:,1);
    TempShallow = DepthShallowA(:,3);
    DatesShallow = DepthShallowA(:,2);
    DatesShallowV = datevec(DatesShallow);
    DatesShallowDay = DatesShallowV(:,1:3);
    DatesShallowDayM = datenum(DatesShallowDay);
    DatesShallowDayMuni = unique(DatesShallowDayM);
    
    binD = floor(min(DepthShallow)):DepthBinIncr:floor(max(DepthShallow));  % depth bins
    Dvec = nan(length(DepthShallow),1);
    for b = 2:length(binD)
        index = find(DepthShallow < binD(b) & DepthShallow > binD(b-1));
        Dvec(index) = binD(b-1);
    end
    binT = floor(min(TempShallow)):TempBinIncr:floor(max(TempShallow)); % temp bins
    Tvec = nan(length(DepthShallow),1);
    for b = 2:length(binT)
        index = find(TempShallow < binT(b) & TempShallow > binT(b-1));
        Tvec(index) = binT(b-1);
    end
end

clear b index

%% figures

for k = 1:length(DatesShallowDayMuni)
    klog = DatesShallowDayMuni(k) == DatesShallowDayM;
    
    
    figure % do subplots
    hist(TempShallow(klog), binT) % this specifies the bins
    TempHist = hist(TempShallow(klog), binT);
    xlabel('Temperature','FontSize', 14, 'FontWeight', 'bold')
    ylabel('Count','FontSize', 14, 'FontWeight', 'bold')
    title({['Combined Tag Temperature Histograms: Shallow (<' num2str(ShallowCutoff) 'm)']; ['Median Temp: ' num2str(ShallowSumm(4))]}, 'FontSize', 14, 'FontWeight', 'bold');
    % legend(num2str(TagIDsUni))
    FigName = ['TempHistosCombinedShallow' num2str(ShallowCutoff) 'm.pdf'];
    annotate_JS(Mfilename, gcf, FigName)
    % orient landscape % save
    % print('-dpdf', FigName);
end

figure
hist(TempAll, 1:length(TagIDsUni)) % this specifies the bins
TempHist = hist(TempShallow, binT);
xlabel('Temperature','FontSize', 14, 'FontWeight', 'bold')
ylabel('Data Counts','FontSize', 14, 'FontWeight', 'bold')
title('Count Contributions from Each Individual Tag', 'FontSize', 14, 'FontWeight', 'bold');  
FigName = 'HistCountEachTag.pdf';
annotate_JS(Mfilename, gcf, FigName)
% orient landscape % save
% print('-dpdf', FigName);

%% 
  
disp('Completed compareTagsSST_EachTag.m') 
% ===== EOF [compareTagsSST_EachTag.m] ======  
