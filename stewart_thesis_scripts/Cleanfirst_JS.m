function [M,Headers] = Cleanfirst_JS(filename)

% CLEANFIRST_JS.m Creates a .mat file with archival tag data. 
% reads an archival .csv file made with Wildlife Computer's Hex Decoder 
% (date set to mm/ddd/yyyy). Can be Very Slow for long deployments.
% 
% changes:
% will import date in format from 2009+ (time series data output), eg:
% "Date","Time","Depth","External Temperature","Light Level"
% 05-Dec-2009,02:59:00,3.0,16.75,157.25
%
% Sample Call: 
%	[M,Headers] = Cleanfirst_JS(filename);
%  OR 
%   [M,Headers] = Cleanfirst_JS; % use if you want to use a browser to find the
%   .csv file
%   
% Inputs: 
%  filename -> file path: 'C:\Data\file.csv'
% Outputs: 
%  M -> is an matlab format of the file
%   the funtion exports a .csv file and a .mat file
%
% Outside Functions Called: 
%   keep.m (MATLAB Central)optional
%  
% AUTHOR    : A. Booth ashley@boothfamily.com, J. Stewart jules32@gmail.com 
% DATE      : 21-Sep-2011 13:39:44, based on A. Booth's CleanFirst_AB.m
% Revision  : 1.00  
% DEVELOPED : 7.4.0.287 (R2007a) Mac 
% FILENAME  : Cleanfirst_JS.m 

%% Import csv file
if nargin == 0
    [csvFile,dir1] = uigetfile('*.csv', 'Select an archival data .csv file made by Wildlife Computer''s Hex Decoder'); % opens a directory browser for user to select file they want
    filename = [dir1 csvFile];
else 
    [pathstr, name, ext, versn] = fileparts(filename);
    csvFile = strcat(name,ext);
end

fid=fopen(filename); %open the file object
tline = fgetl(fid); %grab the first (header) line
fclose(fid); % close file object

%reformat the header line to turn in to a cell array
tline = strrep(tline,'"',''); %if headers have quotes in them, remove them
tline = strrep(tline,' ',''); %if headers have spaces in them, remove them
Tline = ['''' regexprep(tline,',',''';''') '''']; %reformat 
header = eval(['{' Tline '}'])'; %make into string array

ColNum = length(header);
ColNumNum = length(header)-2;

% build format to import data with two columns of strings and a variable
% number of numerical clumns after that
fmat = '%s %s';
for C = 1:ColNumNum
    fmat = strcat(fmat,' %n');
end

%import column data
eval(['[' tline '] = textread(filename, fmat,''headerlines'', 1, ''delimiter'','','');']) 

%set up data for output
numdata = dlmread(filename,',',1,2);

%% Check for bad data

disp('Checking for bad data... ')
TimeL = cellfun(@(x) length(x),Time); %slow
baddex = find(TimeL~=8);
if ~isempty(baddex) %delete bad data rows
    disp('--->> found bad time data')
end

%% Reformat date

if length(Date{1})==11
    year = input(['Type 4-digit deployment year for file ' csvFile ':  ']);
    disp('Reformating dates (may take a few minutes for long files)...')
    tic
    DateTemp = datenum(Date);
    DateV = datevec(DateTemp);
    TimeM = datenum(Time);
    TimeV = datevec(TimeM);
    
    DateM = datenum([DateV(:,1:3) TimeV(:,4:6)]);
    toc
    
% elseif length(Date{1})<11
%     datefmt = input('Type date format in apostrophes (''mm/dd/yyyy''):  ');
%     disp('Reformating dates (may take a few minutes)...')
% %     Time2 = datenum(Time,'HH:MM:SS');
%     DateM =  datenum(Date,datefmt) + (datenum(Time,'HH:MM:SS')-floor(datenum(Time(1),'HH:MM:SS')));
% 
%     
%     [Y, MM, D, H, MN, S] = datevec(DateM);
% else
%     datefmt = input('Type date format in apostrophes (''mm/dd/yyyy HH:MM:SS''):  ');
%     disp('Reformating dates (may take a few minutes)...')
%     DateM = datenum(Date,datefmt);
%     [Y, MM, D, H, MN, S] = datevec(DateM);
else 
    disp('error--check Cleanfirst_JS.m and the .csv file')
end

Headers2 = {header{:},'MATLABDate(UTC)','Year','Month','Day','Hour','Min','Sec'};

%% Prepare for export:
M = [numdata,DateM,[DateV(:,1:3) TimeV(:,4:6)]];
Headers = Headers2(3:end);
%Export to .mat
newfilenameA = strrep(filename,'.csv','');
newfilename = [newfilenameA '.mat'];
save(newfilename,'M','Headers')
disp(['Saved variables to ' newfilename])


disp('Completed Cleanfirst_JS.m') 
% ===== EOF [Cleanfirst_JS.m] ======  
