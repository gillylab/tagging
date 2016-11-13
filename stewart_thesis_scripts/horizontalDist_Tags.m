% horizontalDist_Tags.m
%
% Called From:

% Description:
%
% Outside Functions Called:
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 20-Apr-2011 10:38:48
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : horizontalDist_Tags.m

Mfilename = mfilename;

%% Horizontal distances

distCAorMx = input('distance calculation for CA (1),  GOC (2) of Magda (3)? **for CA tags, make sure all are selected');

if distCAorMx == 1 % CA tags
    DeployDays = [
        2.7    % 83046_08
        12.7   % 83046_09
        5.7    % 83052
        3      % 64004
        17.6]; % 83051
    LocTagDeploy = LocTagDeployCA;
    LocTagPopUp = LocTagPopUpCA;
    
elseif distCAorMx == 2 % GOC tags
    DeployDays = [
        6  %1 8891: datenum('Oct-20-2001')-datenum('Oct-14-2001')
        8  %2 52868: datenum('Nov-02-2004')-datenum('Oct-25-2004')
        8  %3 52869: datenum('Nov-02-2004')-datenum('Oct-25-2004')
        5  %4 52865: datenum('Oct-31-2004')-datenum('Oct-26-2004')
        3  %5 62006: datenum('Nov-14-2005')-datenum('Nov-11-2005')
        4  %6 62007: datenum('Nov-15-2005')-datenum('Nov-11-2005')
        4];%7 62009: datenum('Nov-15-2005')-datenum('Nov-11-2005')
    
    LocTagDeployGOC = [ % GOC tag locations
        -112.21, 27.34; % 8891
        -112.22, 27.34; % 52868
        -112.22, 27.34; % 52869
        -112.24, 27.35; % 52865
        -112.17, 27.30; % 62006
        -112.17, 27.29; % 62007
        -112.17, 27.29]; % 62009
    
    LocTagPopUpGOC = [ % GOC tag locations
        -112.04, 27.33;
        -110.78, 26.29;
        -111.01, 27.43;
        -109.87, 26.22;
        -111.71, 27.90;
        -112.34, 28.16;
        -111.92, 27.50];
    
    LocTagDeploy = LocTagDeployGOC;
    LocTagPopUp = LocTagPopUpGOC;
    
elseif distCAorMx == 3
    DeployDays = [
        9   %1 52910: datenum('Jun-19-2005')-datenum('Jun-10-2005')
        13  %2 54560: datenum('JUn-23-2005')-datenum('Jun-10-2005')
        23];%3 62009: datenum('Jul-03-2005')-datenum('Jun-10-2005')
    
    LocTagDeployMagda = [ % GOC tag locations
        -112.01, 24.45;
        -112.01, 24.45;
        -112.01, 24.45];
    
    LocTagPopUpMagda = [ % GOC tag locations
        -111.80, 22.99;
        -112.31, 24.44;
        -112.24, 23.81];
    
    LocTagDeploy = LocTagDeployMagda;
    LocTagPopUp = LocTagPopUpMagda;
end


DeployHours = DeployDays*24;

Dist_km = [];
Rate_km_hr = [];
Rate_km_da = [];
for i = 1:size(LocTagDeploy,1)
    [dist_deg,az] = distance(LocTagDeploy(i,2), LocTagDeploy(i,1), LocTagPopUp(i,2), LocTagPopUp(i,1));
    dist_km = deg2km(dist_deg);
    rate_km_hr = dist_km/DeployHours(i);
    rate_km_da = rate_km_hr*24;
    
    Dist_km = [Dist_km; dist_km];
    Rate_km_hr = [Rate_km_hr; rate_km_hr];
    Rate_km_da = [Rate_km_da; rate_km_da];
end

Rate_km_hrAll = [Rate_km_hrAll; Rate_km_hr]; % run this for all 3 groups.
if length(Rate_km_hrAll) == 15
    summary(Rate_km_hrAll)
    plot(1:length(Rate_km_hrAll),Rate_km_hrAll,'.')
    hist(Rate_km_hrAll)
end


%%

disp('Completed horizontalDist_Tags.m')
% ===== EOF [horizontalDist_Tags.m] ======
