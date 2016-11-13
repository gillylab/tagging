% compareTagsMBARI.m
%
% Called From:

% Description:
%
% Outside Functions Called:
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 28-Mar-2012 21:29:01
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : compareTagsMBARI.m

Mfilename = mfilename;
DirDepth = '/Users/juliastewart/Dropbox/Thesis/Dg_Abundance/MBARI_Dg/_figures/_figsThesisCh3/_DepthAnalyses/';

%%

matfilefun = @(x) char(regexp(x, '.+HistoDaily.mat','match')); %created in processSeriesData.m

% read in histogram info
matfilename = cellfun(matfilefun,{matfiles(:).name}, 'UniformOutput',0); %find files of interest within cell array
matfiledex = cellfun(@isempty,matfilename);
matfilename = char(matfilename(~matfiledex));
matfilename = cellstr(matfilename);

cc = ['r','b','k','m','c','g','r','b','k','m','c','g','r','b','k','m','c','g',]; % add seomthing at the beginning of CCS-4 to get the color right
ss = ['-kx'; '-ko'; '-k+'; '-ks'; '-kd'];
for i = 1%:length(matfilename)
    filetemp = ['load ' matfilename{i}];
    eval(filetemp);
    figure; hold on
    for k=2:size(DepthDperc_daily,2)-1 % depending on the tag, don't include the first day (i=1,4) or last day (i=6)
        H = plot(dd,DepthDperc_daily(:,k), 'o-','color', cc(k), 'linewidth', 1.5);
        %         H = plot(dd,DepthDperc_daily(:,k), ss(k,:), 'linewidth', 1.5, 'markersize', 10);
        %         legend(H, num2str(k))
    end
    title(matfilename{i})
    set(gca, 'XLim', [0 1000])
    ylabel('Percent of Time','FontSize', 14, 'FontWeight', 'bold')
    xlabel('Depth (m) ', 'FontSize', 14, 'FontWeight', 'bold')
    set(gca, 'fontsize', 12, 'fontweight', 'bold')
    FigName = ['Tag' num2str(i) 'DailyDepthHistPerc.pdf'];
    annotate_JS(Mfilename, gcf, FigName)
    orient landscape % save
    print('-dpdf', [DirDepth FigName]);
    switch i
        case 1
            Tag1depthDperc_daily = DepthDperc_daily;
            Tag1depthbins = dd;
              Tag1maxdepth = [283.626444816602,356.473490771954,409.682442100998,1447.93085291863;];%  % Tag1maxdepth = CA_5tagMaxTempDepth(1:4,6)
        case 2
            Tag2depthDperc_daily = DepthDperc_daily;
            Tag2depthbins = dd;
            Tag2maxdepth = CCSmaxDtemp;
        case 3
            Tag3depthDperc_daily = DepthDperc_daily;
            Tag3depthbins = dd;
        case 4
            Tag4depthDperc_daily = DepthDperc_daily;
            Tag4depthbins = dd;
        case 5
            Tag5depthDperc_daily = DepthDperc_daily;
            Tag5depthbins = dd;
        case 6
            Tag6depthDperc_daily = DepthDperc_daily;
            Tag6depthbins = dd;
    end

end

%% 

% save -mat '/Users/juliastewart/Dropbox/Thesis/Dg_Abundance/MBARI_Dg/CA_4tagsDayHist.mat'...
%     CA_4tagDayBins CA_4tagDayPerc CA_4tagMaxTempDepth...
%     Tag1depthDperc_daily  Tag1depthbins Tag1maxdepth...
%     Tag2depthDperc_daily  Tag2depthbins Tag2maxdepth...
%     Tag3depthDperc_daily  Tag3depthbins...
%     Tag4depthDperc_daily  Tag4depthbins...
%     Tag5depthDperc_daily  Tag5depthbins...
%     Tag6depthDperc_daily  Tag6depthbins


%%

disp('Completed compareTagsMBARI.m')
% ===== EOF [compareTagsMBARI.m] ======
