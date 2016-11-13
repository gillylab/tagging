% plotContoursTime.m
%  
% Outputs: 
%  
% Outside Functions Called: 
%  
% AUTHOR    : J. Stewart jules32@gmail.com 
% DATE      : 18-Mar-2011 17:43:18  
% Revision  : 1.00  
% DEVELOPED : 7.9.0.529 (R2009b) OSX 
% FILENAME  : plotContoursTime.m 

Mfilename = mfilename; 

%%

% depth bins
binsizeD = 5;
dbar = binsizeD:binsizeD:max(DepthTag);
Depth_binnedN = hist(DepthN,dbar); %count within each bin
Depth_binnedD = hist(DepthD,dbar); 
percDN = ((Depth_binnedN/sum(Depth_binnedN))*100)';
percDD = ((Depth_binnedD/sum(Depth_binnedD))*100)';

%% Setup for DAILY figures below: Calc percent time at depth/temp

% setup for daily matrices of histogram output
DepthNbinned_daily = [];
DepthDbinned_daily = [];
DepthNperc_daily = [];
DepthDperc_daily = [];

DeployHours = [];
DepthTagMax = []; % Also get max depth for each day for Hinke model

for i = 1:length(DEPLOYINDEXuni)
    dlog = DEPLOYINDEX == DEPLOYINDEXuni(i);
    
    % daily depth bins
    DepthNbinnedi= hist(DepthN(dlog),dbar)'; %count within each bin    
    DepthDbinnedi= hist(DepthD(dlog),dbar)'; 
    DepthNperci = ((DepthNbinnedi/sum(DepthNbinnedi))*100);
    DepthDperci = ((DepthDbinnedi/sum(DepthDbinnedi))*100);
      
    if exist('pttdata', 'var')
        if i == 1 && pttdata(1) == 83051 % Don't understand why this is giving a problem
            DepthDperci = DepthDperci';
            
        end
    end
    
    %for troubleshooting
    if size(DepthNperci,2) > size(DepthNperci,1) % if it's a row not column vector
        DepthNperci = DepthNperci';
    end
    if size(DepthDperci,2) > size(DepthDperci,1) 
        DepthDperci = DepthDperci';
    end
    
    DepthNbinned_daily = [DepthNbinned_daily DepthNbinnedi]; 
    DepthDbinned_daily = [DepthDbinned_daily DepthDbinnedi];
    DepthNperc_daily = [DepthNperc_daily DepthNperci];
    DepthDperc_daily = [DepthDperc_daily DepthDperci];
  
end

%%
DEPLOYINDEXpoint5 = DEPLOYINDEX(1):0.5:DEPLOYINDEX(end)+0.5;
dp5 = datenum(yr,mon,DEPLOYINDEXpoint5,0,0,0);

xx = [];
for i = 1:length(dp5)
    x = repmat(dp5(i),length(dbar),1);
    xx = [xx x];
end

yy = repmat(dbar,length(DEPLOYINDEXpoint5),1)';

zz = [];
for i = 1:length(DEPLOYINDEXuni)
    zz = [zz DepthNperc_daily(:,i)/100 DepthDperc_daily(:,i)/100];
end
zz(isnan(zz)) = 0;


%% plot

figure
pcolor(xx,yy,zz)
shading flat  %shading interp to make it smooth, but this can clip somewhat, so used flat to troubleshoot
set(gca,'YDir','reverse', 'YLim', [0 1500]);
colormap(flipud(winter))
colorbar
datetick('x', 'keeplimits')
title([tagName ' time-at-depth, ' num2str(binsizeD) 'm bins '],'FontSize',20, 'FontWeight', 'bold')
ylabel('Depth (m)', 'FontSize',16, 'FontWeight', 'bold')
xlabel('Date (PST)  ', 'FontSize',16, 'FontWeight', 'bold')
set(gca, 'fontsize', 12, 'fontweight', 'bold')
FigName = [tagName '_Profile_ContourTimeWinter'];
annotate_JS(Mfilename, gcf, FigName);
% orient landscape % save
% print('-dpdf', [FigName '.pdf']);
% print('-dpng', [FigName '.png']);
% print('-djpeg', [FigName '.jpg']);

figure
pcolor(xx,yy,zz)
shading flat  %shading interp to make it smooth, but this can clip somewhat, so used flat to troubleshoot
set(gca,'YDir','reverse', 'YLim', [0 1500]);
colormap(flipud(bone))
colorbar
datetick('x', 'keeplimits')
title([tagName ' time-at-depth, ' num2str(binsizeD) 'm bins '],'FontSize',20, 'FontWeight', 'bold')
ylabel('Depth (m)', 'FontSize',16, 'FontWeight', 'bold')
xlabel('Date (PST)  ', 'FontSize',16, 'FontWeight', 'bold')
set(gca, 'fontsize', 12, 'fontweight', 'bold')
FigName = [tagName '_Profile_ContourTimeBone'];
% annotate_JS(Mfilename, gcf, FigName);
% orient landscape % save
% print('-dpdf', [FigName '.pdf']);
% print('-dpng', [FigName '.png']);
% print('-djpeg', [FigName '.jpg']);

%% 
  
disp('Completed plotContoursTime.m') 
% ===== EOF [plotContoursTime.m] ======  
