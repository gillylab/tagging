% import2010CTD.m
%  
% Outputs: 
%  
% Outside Functions Called: 
%       convertDOunits.m % A. Booth  
%       sw_pres.m called from CSIRO Seawater MATLAB toolkit: http://www.cmar.csiro.au/datacentre/ext_docs/seawater.htm
% 
% AUTHOR    : J. Stewart jules32@gmail.com 
% DATE      : 02-Feb-2011 09:51:22  
% Revision  : 1.00  
% DEVELOPED : 7.9.0.529 (R2009b) OSX 
% FILENAME  : import2010CTD.m 

%% 

cd '/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/_CTDs/_FulmarVIPtrip_2010_11_30';
filename = 'FulmarVIPNov2010.txt'; % this was created using A. Booth's plotCTD.m, in the processingCTDs_Booth_2011_Jan folder. 

fid=fopen(filename);
line = 1;
tline = fgetl(fid);
fclose(fid);
Tline = ['''' regexprep(tline,'\t',''';''') '''']; %reformat headers
Nov2010CTD.headernames = eval(['{' Tline '}']); %make into string array
Nov2010CTD.headernames = Nov2010CTD.headernames';
numcol = length(Nov2010CTD.headernames);

%import column data
fmat = '%n %s %n %n %n %n %n %n %n %n %n';
fid=fopen(filename);
Nov2010CTD.data = textscan(fid, fmat,'HeaderLines', line, 'Delimiter','\t');
fclose(fid);

%% cleanup

Nov2010date = datenum(Nov2010CTD.data{1,2});

% setup for convertDOunits by A. Booth
DO = Nov2010CTD.data{1,8}; % perctsat
% DO = Nov2010CTD.data{1,9}; % umolkg
Units = 'perctsat';
Temp_C = Nov2010CTD.data{1,11};
Sal_psu = Nov2010CTD.data{1,10};
Depth_m = Nov2010CTD.data{1,6};
Lat = Nov2010CTD.data{1,4};
Press_db = sw_pres(Depth_m,Lat);
[DO_umolkg,DO_umolL,DO_mgL,DO_mlL,DO_perctsat] = convertDOunits(DO,Units,Temp_C,Sal_psu,Press_db); % A. Booth

matrixNov2010CTD = [Nov2010CTD.data{1,1} Nov2010date Nov2010CTD.data{1,4} Nov2010CTD.data{1,3} ...
    Nov2010CTD.data{1,6} Nov2010CTD.data{1,11} DO_mlL Nov2010CTD.data{1,10} Nov2010CTD.data{1,8}];
matrixNov2010CTDheaders = {'divenum' 'dateCTD' 'lat' 'lon' 'depth' 'temper' 'oxyg_umolKg' 'salin' 'oxyg_perc_sat'};

%% save as mat file

keep matrixNov2010CTD matrixNov2010CTDheaders 

cd('/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/_CTDs/_FulmarVIPtrip_2010_11_30');
save -mat Nov2010CTD.mat

%%
  
disp('Completed import2010CTD.m') 
% ===== EOF [import2010CTD.m] ======  
