% plotSeriesDataPublication.m
%
% Called From: plotSeriesData.m

% Description:
%
% Outside Functions Called:
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 29-Nov-2011 12:16:20
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : plotSeriesDataPublication.m

Mfilename = mfilename;
o2 = 1;

%% import oxygen data. Originally made in b_ROVCTD_stats.m and called by calcOMZtimes.m
if o2
    if TagLocation == 1
        %         cd('/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/_CTDs');
        %         load -mat ROVROVfalltib_depth-m_PolyFit_4_12casts_NovDec.mat % All: Casts in November-December 2009
        %
        %         OxyFitN_mlL = polyval(P,DepthN); % F = fit
        %         OxyFitD_mlL = polyval(P,DepthD); % F = fit
        %
        %         % convert
        %         [OxyFitN_umolkg,DO_umolL,DO_mgL,DO_mlL,DO_perctsat] = convertDOunits(OxyFitN_mlL,'ml/L',TempN,[],[]); % A. Booth (2010)
        %         [OxyFitD_umolkg,DO_umolL,DO_mgL,DO_mlL,DO_perctsat] = convertDOunits(OxyFitD_mlL,'ml/L',TempD,[],[]); % A. Booth (2010)
        
        OMZ = -484.5;
        OLZ = -256.9;
    else
%                 cd('/Users/juliastewart/Dropbox/Thesis/Dg_Tagging/_CTDs/_CTD_tagging_cnvfiles/_GOC/') % don't use our CTDs: don't go deep enough
%                 load -mat 'GOCCTD_PolyInfoTrix.mat'
%         
%                 ppp = input('Which GOC tag number? (1-6)');
%                 P = Ptrix(:,ppp);
%                 OxyTag = polyval(P(1:5),DepthTag(DepthTag>15));
        
        OMZ = -237.0; % autumn
        %OMZ = -204.7; %spring
        OLZ = -94.4; % autumn
%         OLZ = -182.6; % spring
    end
end

%% plot

figure
set(gcf,'Position',[12   213   842   688])
scatter(Seriesdate,DepthTag*-1,20,TempTag,'filled');
title([tagName ' Seriesdata2'])
ylabel('Depth (m)')
xlabel('Date (PST)')
Hcbar = colorbar('location','SouthOutside', 'LineWidth', 1);
caxis([4 24]);
xlabel(Hcbar,'Temperature (\circC)');
datetick('x',6)

[x,y] = DayNight_boxes(Seriesdate,sr/24,ss/24,get(gca,'Ylim'),0);
set(gcf,'Position',[12   213   842   688])
fill(x,y,[0.9 0.9 0.9], 'EdgeColor','none')
hold on
scatter(timeN,DepthN*-1,20,TempN,'filled');
scatter(timeD,DepthD*-1,20,TempD,'filled');
caxis([4 24])%colorsLoc); % defined in importTagData_csv.m
plot(timeD, repmat(OMZ,length(timeD),1), '--',timeD, repmat(OLZ,length(timeD),1), '--');
hold off
title([tagName ' ' buflab],'FontSize',20, 'FontWeight', 'bold')
ylabel('Depth (m)', 'FontSize',16, 'FontWeight', 'bold')
xlabel('Day of Deployment ', 'FontSize',16, 'FontWeight', 'bold')
set(gca, 'XLim', [min(Seriesdate) max(Seriesdate)], 'YLim', [-800 0])%[max([DepthN; DepthD])*-1 0]) % min(Seriesdate)+2.71])% for deep dive fig GOC-4
%     set(gca, 'XLim', [min(Seriesdate) max(Seriesdate)-2.25], 'YLim',[-800]) %only for GOC-5
set(gca, 'fontsize', 12, 'fontweight', 'bold')
Hcbar = colorbar('location','SouthOutside', 'LineWidth', 1);
xlabel(Hcbar,'Temperature (\circC)', 'FontSize',16, 'FontWeight', 'bold');
set(gca,'XTick', DeployDates+ datenum(0, 0, 0, 8, 0, 0))
datetick('x','mm/dd', 'keeplimits','keepticks')
FigName = [tagName '_TimeSeries_TimeDepth_Temp' buflab 'PubTalk2.pdf'];
annotate_JS(Mfilename, gcf, FigName);
% orient landscape % save
% print('-dpdf', [dir1 FigName]);

%% plot grayscale
figure
set(gcf,'Position',[12   213   842   688])
scatter(Seriesdate,DepthTag*-1,20,TempTag,'filled');
title([tagName ' Seriesdata2'])
ylabel('Depth (m)')
xlabel('Date (PST)')
Hcbar = colorbar('location','SouthOutside', 'LineWidth', 1);
caxis(colorsLoc);
xlabel(Hcbar,'Temperature (\circC)');
datetick('x',6)

[x,y] = DayNight_boxes(Seriesdate,sr/24,ss/24,get(gca,'Ylim'),0);
set(gcf,'Position',[12   213   842   688])
fill(x,y,[0.9 0.9 0.9], 'EdgeColor','none')
hold on
scatter(timeN,DepthN*-1,20,TempN,'filled');
scatter(timeD,DepthD*-1,20,TempD,'filled');
% colormap(gray)
caxis([0 20])%colorsLoc); % defined in importTagData_csv.m
% scatter(timeN(ooN),DepthN(ooN)*-1,60,OxyFitN_umolkg(ooN));
% scatter(timeD(ooD),DepthD(ooD)*-1,60,OxyFitD_umolkg(ooD));
% plot(timeD, repmat(min(DepthD(ooD))*-1,length(timeD),1), '--');
% colormap(flipud(gray))
% colormap(gray)
hold off
title([tagName ' ' buflab],'FontSize',20, 'FontWeight', 'bold')
ylabel('Depth (m)', 'FontSize',16, 'FontWeight', 'bold')
xlabel('Date (PST)', 'FontSize',16, 'FontWeight', 'bold')
set(gca, 'XLim', [min(Seriesdate) max(Seriesdate)])
set(gca, 'fontsize', 12, 'fontweight', 'bold')
Hcbar = colorbar('location','SouthOutside', 'LineWidth', 1);
% caxis([4 20])%colorsLoc); % defined in importTagData_csv.m
xlabel(Hcbar,'Temperature (\circC)', 'FontSize',16, 'FontWeight', 'bold');
set(gca,'XTick', DeployDates+ datenum(0, 0, 0, 8, 0, 0))
datetick('x','mm/dd', 'keeplimits','keepticks')
FigName = [tagName '_TimeSeries_TimeDepth_Temp' buflab 'pub.pdf'];
annotate_JS(Mfilename, gcf, FigName);
%         orient landscape % save
%         print('-dpdf', [dir1 FigName]);


%%

disp('Completed plotSeriesDataPublication.m')
% ===== EOF [plotSeriesDataPublication.m] ======
