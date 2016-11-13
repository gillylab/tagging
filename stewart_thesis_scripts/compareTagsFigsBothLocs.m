% compareTagsFigsBothLocs.m
%
% Called From: compareTags.m

% Description:
%
% Outside Functions Called:
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 11-Nov-2011; 18-Jan-2012: moved OMZ event durations here and
% made this a stand-alone m file (can clear all before you start)
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : compareTagsFigsBothLocs.m

Mfilename = mfilename;
clear H
dirCompareFigsBoth = '/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/CompareTagDir/_CA_GOC/';
cd(dirCompareFigsBoth);

%% Make histograms with both locations together

load -mat 'AllTagsCompareDataCA.mat'; % made in compareTags.m
BPdepthCA = BoxPlotDepthtrix;
BPtempCA = BoxPlotTemptrix;

HdepthCA = DepthHistoTrix;
HtempCA = TempHistoTrix;
HdepthLengthCA = DepthHistotrixLength;
ID_CA = TagIDs;


load -mat 'AllTagsCompareDataGOC.mat';
BPdepthGOC = BoxPlotDepthtrix;
BPtempGOC = BoxPlotTemptrix;
BPdepthGOC(:,5) = BPdepthGOC(:,5)+100; % label tags differently
BPtempGOC(:,5) = BPtempGOC(:,5)+100;

HdepthGOC = DepthHistoTrix;
HtempGOC = TempHistoTrix;
HdepthLengthGOC = DepthHistotrixLength;
ID_GOC = TagIDs+100;

dbar = 0:10:max([BPdepthCA(:,3);BPdepthCA(:,4)]); % use CA because deepest dive
tbar = (floor(min([BPtempCA(:,3);BPtempCA(:,4)])):0.5:max([BPtempGOC(:,3);BPtempGOC(:,4)])); % use widest temp range


for p = 1%:2 % Depth or Temp
    if p == 1 % Depth % this is going to have a day and night component
        Data1 = BPdepthCA;
        Data2 = BPdepthGOC;
        DataLabel = 'Depth-m';
        ydir = 'reverse';
        %         percCA = HdepthCA(:,4); % setup for histos
        %         percGOC = HdepthGOC(:,5);
        Bar = dbar;
        %         Bin = dbar(2) - dbar(1);
        barend = max([HdepthLengthCA HdepthLengthGOC]); % for ca: DepthHistotrixLength(1) aug 1 2011
        %         dbarAll = 0:Bin:Bin*barend;
        DataXLim = 800;%1500;
        DataYLim = 25;%40;
    elseif p == 2 % Temp
        Data1 = BPtempCA;
        Data2 = BPtempGOC;
        DataLabel = 'Temperature-C';
        ydir = 'normal';
        %         percCA = HtempCA(:,4); % setup for histos
        %         percGOC = HtempGOC(:,5);
        Bar = tbar;
        %         Bin = dbar(2) - dbar(1);
        barend = max([HdepthLengthCA HdepthLengthGOC]); % for ca: TempHistotrixLength(1) aug 1 2011
        %         dbarAll = 0:Bin:max(dbar);
        DataXLim = 30;
        DataYLim = 30;
    end
    
    
    Bin = Bar(2) - Bar(1);
    barend = max(Bar);
    dbarAll = 0:Bin:max(Bar);
    
    TagIDstats = []; % for the table below
    TagAllstats = [];
    
    % setup for boxplots
    DataDayCA = Data1(:,4);
    DataNightCA = Data1(:,3);
    DataDayGOC = Data2(:,4);
    DataNightGOC = Data2(:,3);
    DatamaxD = max([DataDayGOC; DataDayCA]);
    DatamaxN = max([DataNightGOC; DataNightCA]);
    
    DataDay = [DataDayCA; DataDayGOC];
    DataNight = [DataNightCA; DataNightGOC];
    DataID = [Data1(:,5); Data2(:,5)];
    %
    figure
    subplot(1,2,1)
    boxplot(DataDay, DataID, 'colors', 'k','symbol', '.k', 'outliersize', 5); % plot 5 tags
    %     boxplot([DataCA; DataCA], [Data(:,5); repmat(7,length(DataCA),1)], 'colors', 'k','symbol', '.k', 'outliersize', 3); % plot individs and all tags
    set(gca, 'YLim', [0-Bin DataXLim], 'YDir', ydir) % since boxplots are still on y axis
    title('Nighttime','FontSize', 14, 'FontWeight', 'bold')
    subplot(1,2,2)
    boxplot(DataNight, DataID, 'colors', 'k', 'symbol', '.k', 'outliersize', 5); % plot 5 tags
    %     boxplot([DataGOC; DataGOC], [Data(:,5); repmat(7,length(DataGOC),1)], 'colors', 'k','symbol', '.k', 'outliersize', 3); % plot individs and all tags
    set(gca, 'YLim', [0-Bin DataXLim],'YDir', ydir)
    title('Daytime','FontSize', 14, 'FontWeight', 'bold')
    
    key = [DataLabel ' All Boxplots CCS-GOC'];
    [h3] = suptitle(key);
    [ax,h(1)] = suplabel('Tag Number ','x'); %('Day of the Year ','x');
    [ax,h(2)] = suplabel(DataLabel, 'y');
    set(h3,'FontSize',24, 'FontWeight', 'bold')
    set(h,'FontSize',20, 'FontWeight', 'bold')
    FigName = ['BoxplotSubplot_AllTags' DataLabel '_CCS-GOC.pdf'];
    annotate_JS(Mfilename, gcf, FigName)
    fixfig
%     orient landscape % save
%     print('-dpdf', FigName);
    
    
    %%%% histograms%%%%%
    % combined histograms for depth and temperature setup
    for q = 1:2
        
        if q == 1 % day
            HistDataCA = DataDayCA;
            HistDataGOC = DataDayGOC;
            dnlabel = 'Day';
        else % night
            HistDataCA = DataNightCA;
            HistDataGOC = DataNightGOC;
            dnlabel = 'Night';
        end
        
        Data_binnedC = hist(HistDataCA,dbarAll);%count within each bin
        Data_binnedG = hist(HistDataGOC,dbarAll);
        percCA = ((Data_binnedC/sum(Data_binnedC))*100)';
        percGOC = ((Data_binnedG/sum(Data_binnedG))*100)';
        
        summdepthGOC = summary(HistDataGOC);
        summdepthCA = summary(HistDataCA);
        
        figure;
        H(2) = bar(dbarAll,percGOC,'BarWidth',1);
        hold on
        H(1) = bar(dbarAll,percCA,'BarWidth',1);
        hold off
        set(H(1),'facecolor','none','EdgeColor','k','LineWidth',1.5) %night
        set(H(2),'facecolor',[0.5020    0.5020    0.5020],'EdgeColor',[0.2 0.2 0.2],'LineWidth',1.5) %day
        set(gca, 'XLim', [0-Bin DataXLim], 'YLim', [0 DataYLim], 'XMinorTick', 'on')
        ylabel('Percent of Time','FontSize', 14, 'FontWeight', 'bold')
        xlabel(DataLabel,'FontSize', 14, 'FontWeight', 'bold')
        legend(H,'CCS','GOC','Location','NorthEast')
        title({[' Combined Tag Histograms: CCS-GOC '  DataLabel ' '  dnlabel ' ' ]; ['GOC mean = ' num2str(summdepthGOC(1))...
            ', sd = ' num2str(nanstd(HistDataGOC)) ', median = ' num2str(summdepthGOC(4))];...
            ['CCS mean = ' num2str(summdepthCA(1)) ', sd = ' num2str(nanstd(HistDataCA))...
            ', median = ' num2str(summdepthCA(4)) ]});
        FigName = ['HistosAllTags_' DataLabel '_CCS-GOC_' dnlabel 'shallow.pdf']; % when run with just forage tags, label it.
        annotate_JS(Mfilename, gcf, FigName)
%                 orient landscape % save
%                 print('-dpdf', FigName);

        DirDefense = '/Users/juliastewart/Dropbox/Stanford/2012things/_thesisSubmission/_keynote/';
        % TO GET THIS TO RUN, change dbar above to 25. 
        figure;
        H(2) = plot(percGOC,dbarAll, 'o-k', 'linewidth', 1);
        hold on
        H(1) = plot(percCA,dbarAll, 'o-b', 'linewidth', 1);
        hold off
        set(gca, 'YLim', [0 800], 'XLim', [0 100], 'XMinorTick', 'on', 'YDir', 'reverse')
        xlabel('Percent of Time','FontSize', 14, 'FontWeight', 'bold')
        ylabel(DataLabel,'FontSize', 14, 'FontWeight', 'bold')
        fixfig
        legend(H,'CCS','GOC','Location','NorthEast')
        title({[' Combined Tag Histograms: CCS-GOC '  DataLabel ' '  dnlabel ' ' ]; ['GOC mean = ' num2str(summdepthGOC(1))...
            ', sd = ' num2str(nanstd(HistDataGOC)) ', median = ' num2str(summdepthGOC(4))];...
            ['CCS mean = ' num2str(summdepthCA(1)) ', sd = ' num2str(nanstd(HistDataCA))...
            ', median = ' num2str(summdepthCA(4)) ]});
        FigName = ['HistosAllTags_' DataLabel '_CCS-GOC_' dnlabel 'LINES.pdf']; % when run with just forage tags, label it.
        annotate_JS(Mfilename, gcf, FigName)
%                 orient landscape % save
%                 print('-dpdf', [DirDefense FigName]);

    end
end

%% OXYGEN 

%%%%%setup for histograms%%%%
load -mat ROVROVfalltib_depth-m_PolyFit_4_12casts_NovDec.mat % for CCS CTD info. made in calcOMZtimesCCS
load -mat 'GOCCTD_PolyInfoTrix.mat' % for GOC CTD info. made in calcOMZtimesGOC
binsizeO = 5;
DataXLim = 300;
DataYLim = 20;

%%%%% setup for event durations%%%%%%
load -mat 'OMZeventsPolynomial_hr_CA.mat' % load CA. made in calcOMZtimesCCS. Also can override and save as OMZeventsPolynomial_hr_CA_ALLDAYS.mat for comparison.
% load -mat 'OMZ60eventsPolynomial_hr_CA.mat' %this is if you want to run with it from 20-60 umol/kg
durH_CA = inOMZ_H;
durM_CA = inOMZ_M;
load -mat 'OMZeventsPolynomial_hr_GOC.mat' % load GOC. made in calcOMZtimesGOC
% load -mat 'OMZ60eventsPolynomial_hr_GOC.mat' % load GOC. made in calcOMZtimesGOC
durH_GOC = inOMZ_H;
durM_GOC = inOMZ_M;

% set either mins or hours
MM = 1; % MM=1 for minutes, MM=0 for hours
if MM
    durEvents_CA = durM_CA;
    durEvents_GOC = durM_GOC;
    LL = 'minutes';
else
    durEvents_CA = durH_CA;
    durEvents_GOC = durH_GOC;
    LL = 'hours';
end

tCA = durEvents_CA(:); % temporary
plotCA = tCA(tCA ~= 0);
tGOC = durEvents_GOC(:); % temporary
plotGOC = tGOC(tGOC ~= 0);

%% plot! 

%%%% histograms%%%%%
for q = 1:2 % plot both temperature and depth
    if q == 1 % daytime
        OxyFitGOC_umolkg = OxyFitD_individsGOC_umolkg(:,1); % GOC
        OxyFitCA_umolkg = OxyFitD_tagsCA_umolkg; % CCS
        dnlabel = 'Day';
    else % nighttime
        OxyFitGOC_umolkg = OxyFitN_individsGOC_umolkg(:,1); % GOC
        OxyFitCA_umolkg = OxyFitN_tagsCA_umolkg; % CCS
        dnlabel = 'Night';
    end
    
    obar = 0:binsizeO:max([OxyFitGOC_umolkg;OxyFitCA_umolkg]);
    Oxy_binnedG = hist(OxyFitGOC_umolkg,obar); %count within each bin
    Oxy_binnedC = hist(OxyFitCA_umolkg,obar); %count within each bin
    percOG = ((Oxy_binnedG/sum(Oxy_binnedG))*100)';
    percOC = ((Oxy_binnedC/sum(Oxy_binnedC))*100)';
    summdepthGOC = summary(OxyFitGOC_umolkg);
    summdepthCA = summary(OxyFitCA_umolkg);
    
    figure;
    H(2) = bar(obar,percOG,'BarWidth',1);
    hold on
    H(1) = bar(obar,percOC,'BarWidth',1);
    H(3) = plot(repmat(20,DataYLim+1, 1),0:DataYLim, '--k', 'linewidth', 2)
    H(4) = plot(repmat(60,DataYLim+1, 1),0:DataYLim, '-.k', 'linewidth', 2);
    hold off
    set(H(1),'facecolor','none','EdgeColor','k','LineWidth',1.5) %cal
    set(H(2),'facecolor',[0.5020    0.5020    0.5020],'EdgeColor',[0.2 0.2 0.2],'LineWidth',1.5) %goc
    set(gca, 'XLim', [0-Bin DataXLim], 'YLim', [0 DataYLim], 'XMinorTick', 'on')
    ylabel('Percent of Time','FontSize', 14, 'FontWeight', 'bold')
    xlabel('Oxygen (umol/kg)', 'FontSize', 14, 'FontWeight', 'bold')
    legend(H,'CCS','GOC','Location','NorthEast')
    title({[' Combined Tag Histograms: CCS-GOC Oxygen: ' dnlabel  ' ']; ['GOC mean = ' num2str(summdepthGOC(1))...
        ', sd = ' num2str(nanstd(OxyFitGOC_umolkg)) ', median = ' num2str(summdepthGOC(4))];...
        ['CCS mean = ' num2str(summdepthCA(1)) ', sd = ' num2str(nanstd(OxyFitCA_umolkg))...
        ', median = ' num2str(summdepthCA(4)) ]});
    FigName = ['HistosAllTags_Oxygen_CCS-GOC_' dnlabel '.pdf']; % By hand saved it as ALLDAYS after hand-running OMZeventsPolynomial_hr_CA_ALLDAYS.mat
    annotate_JS(Mfilename, gcf, FigName)
    %         orient landscape % save
    %         print('-dpdf', FigName);
end


% plot event durations
figure; hold on
[n, xout] = hist(plotGOC, max(plotGOC)/10);
[n2, xout2] = hist(plotCA,xout); % for stacking
N = [n2; n; nan(1,length(n2))]'; % to make the colormap work with stacking
% H = bar(xout, N, 'stacked');
% subplot(2,1,1)
H = plot(xout(n>0), n(n>0), 'ro', xout(n2>0), n2(n2>0), 'ko');
% colormap(flipud(gray))
% set(H(1),'facecolor','none','EdgeColor','k','LineWidth',1.5) %night
% set(H(2),'facecolor',[0.5020    0.5020    0.5020],'EdgeColor',[0.2 0.2 0.2],'LineWidth',1.5) %day
xlabel({['Duration below OMZ depth (' LL ')']},'FontSize', 14, 'FontWeight', 'bold')
ylabel('Number of Events','FontSize', 14, 'FontWeight', 'bold')
legend(H,'CCS','GOC','Location','NorthEast')
title({['Duration of Events in OMZ ']; [num2str(length(plotGOC)) ' GOC events']; [num2str(length(plotCA)) ' CCS events']},...
    'FontSize', 14, 'FontWeight', 'bold');
FigName = ['OMZevents_GOC_CCS.pdf']; % add a 60 if running that
annotate_JS(Mfilename, gcf, FigName)
% orient landscape % save
% print('-dpdf', FigName);


% repeat, just for making figure more pretty in Illustrator:
figure; hold on
% H = bar(xout(2:end), N(2:end,:), 'stacked');
H = plot(xout(n>0), n(n>0), '+k', xout(n2>0), n2(n2>0), 'ko');
colormap(gray)
set(gca, 'YLim', [0 25], 'XMinorTick', 'on', 'YMinorTick', 'on') % 25 or 60
xlabel({['Duration below OMZ depth (' LL ')']},'FontSize', 14, 'FontWeight', 'bold')
ylabel('Number of Events','FontSize', 14, 'FontWeight', 'bold')
FigName = ['OMZevents_GOC_CCS-ai.pdf']; % add a 60 if running that
annotate_JS(Mfilename, gcf, FigName)
% orient landscape % save
% print('-dpdf', FigName);

% a few stats for durations
timeset = 10;
shortGOClog = plotGOC <= timeset;
shortGOC = sum(shortGOClog)/length(plotGOC)*100 
shortCAlog = plotCA <= timeset;
shortCA = sum(shortCAlog)/length(plotCA)*100 

longGOClog = plotGOC >= timeset*12; % 10*12= 2 hours
longGOC = sum(longGOClog)/length(plotGOC)*100 
longCAlog = plotCA >= timeset*12;
longCA = sum(longCAlog)/length(plotCA)*100 

%%

disp('Completed compareTagsFigsBothLocs.m')
% ===== EOF [compareTagsFigsBothLocs.m] ======
