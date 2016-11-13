% subsampleArchivalTag.m
%
% Called From: importTagData_csv.m

% Description: subsample archival data from 1 second resolution to 75
% second resolution. Does not average, just pulls the value from every 75th
% one.
%
% Outside Functions Called:
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 27-Jul-2011 14:18:46
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : subsampleArchivalTag.m


%% setup

lengthM = length(DepthTag)/sampleint;
DateRaw = nan(floor(lengthM),1);
DepthRaw = nan(floor(lengthM),1);
TempRaw = nan(floor(lengthM),1);
for i = 1:floor(lengthM)
    DateRaw(i) = Seriesdate(sampleint*i);
    DepthRaw(i) = DepthTag(sampleint*i);
    TempRaw(i) = TempTag(sampleint*i);
end
Seriesdate = DateRaw;
DepthTag = DepthRaw;
TempTag = TempRaw;

clear sampleint lengthM DateRaw DepthRaw TempRaw

% To save in importTagData_csv. Removed this from processSeriesData at the
% bottom
% if TagRecovered
%     sampleint = 75;
%     lengthM = length(DepthN)/sampleint;
%     DateRaw = nan(floor(lengthM),1);
%     DepthRaw = nan(floor(lengthM),1);
%     TempRaw = nan(floor(lengthM),1);
%     BoxplotrixRaw = nan(floor(lengthM),size(Boxplotrix,2));
%
%     for i = 1:floor(lengthM)
%         DateRaw(i) = Seriesdate(sampleint*i); %%%%%%CHANGE
%         DepthRaw(i) = DepthTag(sampleint*i);
%         TempRaw(i) = TempTag(sampleint*i);
%         BoxplotrixRaw(i,:) = Boxplotrix(sampleint*i,:);
%     end
%     Boxplotrix = BoxplotrixRaw;
% else
%
% end
%
% TagRaw = repmat(tagNumD, length(DateRaw), 1);
% DateRawV = datevec(DateRaw);
%
% clear i sampleint lengthM
%%

disp('Completed subsampleArchivalTag.m')
% ===== EOF [subsampleArchivalTag.m] ======
