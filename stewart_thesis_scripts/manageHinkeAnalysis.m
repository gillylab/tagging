
% load all data from tag [a lot of irrelevant things here]
% load -mat Data_83051.mat % example tag: our Mexico swimmer

habitatID_eachTag % polyfit and PHUI maps. Replaces habitatID. 
% Will prompt user to identify which tag to read in 
%(83051 is the Mexico tag) and then will load data from: 
%         get_stewart
%         tagDailyHistos_Shallow
% the code that created these data are below. habitatID_eachTag right now
% is set up to run for 50m, but can be changed by changing the file to load
% in. 

%% these scripts generate the .mat files that are then run by habitatID_eachTag

get_stewart % set up to use extracto_3D_bdap for Mexico tag


 