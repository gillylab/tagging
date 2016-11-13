% manageTagData.m
%  
% manages tag data/calls other scripts for PAT tags analyzed in:
%	Stewart JS. 2012 Humboldt squid in the northern California Current System. PhD dissertation, Stanford University.
%   Stewart JS, Hazen EL, Foley DG, Bograd SJ, Gilly WF. Marine predator migration during range expansion. Marine Ecology Progress Series.
%   Stewart JS, Field JC, Markaida U, Gilly WF. Deep Sea Research Part II.
%   Behavioral ecology of jumbo squid in relation to oxygen minimum zones.
%  
% AUTHOR    : J. Stewart jules32@gmail.com 
% DATE      : 09-Aug-2010 21:10:37  
% Revision  : 1.00  
% DEVELOPED : 7.9.0.529 (R2009b) OSX 
% FILENAME  : manageTagData.m

%strfind

%% Individual tag analysis:
    %for each individual tag, imports the data, makes and saves figures as well as .mat files in tag's folder

importTagData_csv
        % User input required for tagging location (CA or GOC) and data type (timeseries or PDT)
        % runs processSeriesData and plotSeriesData
        % or runs processPDTData and plotPDTData
        % runs plotROVCTD_Tags, which plots CTD dat over tag data, but also
        % can be run separately. 
        % also calls saveTagDataAsMat, which allows data to be processed by A. Booth's scripts

%% Run Movement Model

habitatID_eachTag 
        % loads files made in get_stewart.m and tagDailyHistosShallow.m (from processSeriesData).
% this runs the movement model but gives user the choice. Calls manageTagData_CTD, which looks at CTD casts in the area of tagging

manageTagData_CTD

%% compare

% compare
compareTags % this processes the import_Series_csv_JS tag *minmaxTD.mat files. This script also calls plotROVCTD_Tags

%%
  
disp('Completed manageTagData.m') 
% ===== EOF [manageTagData.m] ======  
