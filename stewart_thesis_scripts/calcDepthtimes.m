% calcDepthtimes.m
%
% Called From: compareTags.m

% Description: looks at porportion of a 24-hour day that the squid spends
% above 10m, and around OMZ depths...
% really mirrors coding in calcOMZtimes, but this is for depth swaths.
%
% Outside Functions Called:
% prepCTD.m J. Stewart
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 03-Oct-2011 16:04:32
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : calcDepthtimes.m

Mfilename = mfilename;
clear H % so legend won't crash

just1window = 1;
plotthis = 1;

%%

cd(dirBoth)

%%

if ~just1window % do multiple subplots
    
    DepthWindow = [0:25:975; 25:25:1000]'; %give a depth window
    
    for I = 1:size(DepthWindow,1)
        
        Dprop = zeros(DL,length(TagIDs)); % Proportion of Depths above DepthWindow
        Tmed = zeros(DL,length(TagIDs)); % median temperatures corresponding to depths
        Tmean = zeros(DL,length(TagIDs)); % mean
        Tsd = zeros(DL,length(TagIDs)); % for error bars
        Tvar = zeros(DL,length(TagIDs)); % for error bars
        
        for i = 1:length(TagIDs)
            ilog = TagIDs(i) == Rawtrix(:,6);
            DDunique = unique(Rawtrix(ilog,5)); % DeployDays for TagIDs(i)
            for j = 1:length(DDunique)
                jlog = DDunique(j) == Rawtrix(:,5);
                ij = ilog + jlog;
                ijlog = ij == 2;
                Dij = Rawtrix(ijlog,2); % depth
                Tij = Rawtrix(ijlog,3); % temp
                dlog = Dij <= DepthWindow(I,2) & Dij > DepthWindow(I,1); % from compareTags.m: Tag ID that won't change.
                Tijd = Tij(dlog);
                
                % percent in depth window
                depthCount = sum(dlog); % counts
                depthS = depthCount*75; % 75-sec sampling interval (subsampled if archival)
                depthM = depthS/60;
                depthH = depthS/3600; % daytime, in hours
                
                dayij = sum(ijlog); % counts of the deployment (will be just 3 days if that's what the user specified in compareTags.m
                dayijS = dayij*75; % 75-sec sampling interval (subsampled if archival)
                dayijM = dayijS/60;
                dayijH = dayijS/3600;
                
                Dprop(j,i) = depthM/dayijM;
                Tmed(j,i) = nanmedian(Tijd);
                Tmean(j,i) = nanmean(Tijd);
                Tsd(j,i) = std(Tijd);
                Tvar(j,i) = var(Tijd);
                
            end
        end
        
        % reformat
        Dt = Dprop(:); % temporary
        Dmedt = Tmed(:);
        Dmeant = Tmean(:);
        Dsdt = Tsd(:);
        Dvart = Tvar(:);
        
        PropDepths = Dt(Dt > 0);
        PropMedianTemps = Dmedt(Dt > 0);
        PropMeanTemps = Dmeant(Dt > 0);
        PropTempSD = Dsdt(Dt > 0);
        PropTempVar = Dvart(Dt > 0);
        
        % save
%                 filenameSave = ['DepthProportion_' locID '_' num2str(DepthWindow(I,1)) '-' num2str(DepthWindow(I,2)) 'm.mat'];
%                 eval (['save -mat ' filenameSave ' PropDepths', ' PropMedianTemps', ' PropMeanTemps', ' PropTempSD', ' PropTempVar'])
    end
    
    %% load them all and plot
    
    Count = 0;
    numberofpanels = 6;
    numcols = 3;
    if plotthis
        for I = 1:size(DepthWindow,1)
            
            filenameLoad = ['DepthProportion_CA_' num2str(DepthWindow(I,1)) '-' num2str(DepthWindow(I,2)) 'm.mat']; % load CA
            eval (['load -mat ' filenameLoad])
            propD_CA = PropDepths;
            Tmedian_CA = PropMedianTemps;
            Tmean_CA = PropMeanTemps;
            Tsd_CA = PropTempSD;
            Tvar_CA = PropTempVar;
            
            filenameLoad = ['DepthProportion_GOC_' num2str(DepthWindow(I,1)) '-' num2str(DepthWindow(I,2)) 'm.mat']; % load GOC
            eval (['load -mat ' filenameLoad])
            propD_GOC = PropDepths;
            Tmedian_GOC = PropMedianTemps;
            Tmean_GOC = PropMeanTemps;
            Tsd_GOC = PropTempSD;
            Tvar_GOC = PropTempVar;
            
            Count = Count + 1;
            fignum = floor((Count - 1)./numberofpanels)+1; % set the number for 6-panel figures
            panelnum = rem(Count - 1, numberofpanels)+1; % rem takes the remainder after division
            if panelnum == 1
                figure(fignum); clf;
            end
            subplot(ceil(numberofpanels/numcols), numcols, panelnum); hold on;
            
            errorbar_x(Tmean_CA, propD_CA, Tsd_CA, '.')
            H(1) = plot(Tmean_CA, propD_CA, '.','markersize', 15);
            errorbar_x(Tmean_GOC, propD_GOC, Tsd_GOC, '.')
            H(2) = plot(Tmean_GOC, propD_GOC, '.r', 'markersize', 15);
            title({[num2str(DepthWindow(I,1)) '-' num2str(DepthWindow(I,2)) 'm ']}, 'FontSize', 14, 'FontWeight', 'bold');
            set(gca, 'YLim', [0 1], 'XLim', [2 30]);
            set(gca, 'fontsize', 12, 'fontweight', 'bold')
            h1 = gca;
            legend(H,{['CCS, n=' num2str(length(Tmean_CA))],['GOC, n=' num2str(length(Tmean_GOC))]},'Location','NorthEast')
            if panelnum == numberofpanels
                key = {['Proportion of day: depth and temps. ' num2str(fignum)]};
                [h3] = suptitle(key);
                [ax,h(1)] = suplabel(['Mean temperature and SD (°C) ' ],'x');
                [ax,h(2)] = suplabel(['Proportion of 24-hour day '],'y');
                set(h3,'FontSize',24, 'FontWeight', 'bold')
                set(h,'FontSize',20, 'FontWeight', 'bold')
                legend(H,{['CCS, n=' num2str(length(Tmean_CA))],['GOC, n=' num2str(length(Tmean_GOC))]},'Location','NorthEast')
                FigName = ['DepthProp_GOC_CCS_' num2str(DepthWindow(I,2)) 'm.pdf'];
                annotate_JS(Mfilename, gcf, FigName)
%                 orient landscape % save
%                 print('-dpdf', [dirBoth FigName]);
            end
            
        end
        if panelnum < numberofpanels % catch for when there aren't 6 panels on a figure
            key = {['Proportion of day: depth and temps. ' num2str(fignum)]};
            [h3] = suptitle(key);
            [ax,h(1)] = suplabel(['Mean temperature and SD (°C) ' ],'x');
            [ax,h(2)] = suplabel(['Proportion of 24-hour day '],'y');
            set(h3,'FontSize',24, 'FontWeight', 'bold')
            set(h,'FontSize',20, 'FontWeight', 'bold')
            FigName = ['DepthProp_GOC_CCS_' num2str(DepthWindow(I,2)) 'm.pdf'];
            annotate_JS(Mfilename, gcf, FigName)
%             orient landscape % save
%             print('-dpdf', [dirBoth FigName]);
        end
    end
    
else % run for just one depth window
    
    DepthWindow = [0 25]; %give a depth window
    
    Dprop = nan(22,length(TagIDs)); % Proportion of Depths above DepthWindow
    Tmed = nan(22,length(TagIDs)); % median temperatures corresponding to depths
    Tmean = nan(22,length(TagIDs)); % mean
    Tsd = nan(22,length(TagIDs)); % for error bars
    Tvar = nan(22,length(TagIDs)); % for error bars
    
    for i = 1:length(TagIDs)
        ilog = TagIDs(i) == Rawtrix(:,6);
        DDunique = unique(Rawtrix(ilog,5)); % DeployDays for TagIDs(i)
        for j = 1:length(DDunique)
            jlog = DDunique(j) == Rawtrix(:,5);
            ij = ilog + jlog;
            ijlog = ij == 2;
            Dij = Rawtrix(ijlog,2); % depth
            Tij = Rawtrix(ijlog,3); % temp
            dlog = Dij <= DepthWindow(2) & Dij > DepthWindow(1); % from compareTags.m: Tag ID that won't change.
            Tijd = Tij(dlog);
            
            %         figure; hold on
            %         plot(Tij, Dij, '.', 'markersize', 15);
            %         title(num2str(TagNameN(i)))
            %         set(gca,'YDir', 'reverse');
            %         grid on
            
            % percent in depth window
            depthCount = sum(dlog); % counts
            depthS = depthCount*75; % 75-sec sampling interval (subsampled if archival)
            depthM = depthS/60;
            depthH = depthS/3600; % daytime, in hours
            
            dayij = sum(ijlog); % counts of the deployment (will be just 3 days if that's what the user specified in compareTags.m
            dayijS = dayij*75; % 75-sec sampling interval (subsampled if archival)
            dayijM = dayijS/60;
            dayijH = dayijS/3600;
            
            Dprop(j,i) = depthM/dayijM;
            Tmed(j,i) = nanmedian(Tijd);
            Tmean(j,i) = nanmean(Tijd);
            Tsd(j,i) = std(Tijd);
            Tvar(j,i) = var(Tijd);
            
        end
    end
    
    % reformat
    Dt = Dprop(:); % temporary
    Dmedt = Tmed(:);
    Dmeant = Tmean(:);
    Dsdt = Tsd(:);
    Dvart = Tvar(:);
    
    PropDepths = Dt(Dt > 0); %(Dt > 0); % we want all plotted each day--0's are important (although they won't be plotted since no associated temp)
    PropMedianTemps = Dmedt(Dt > 0);
    PropMeanTemps = Dmeant(Dt > 0);
    PropTempSD = Dsdt(Dt > 0);
    PropTempVar = Dvart(Dt > 0);
    
    % save
%         filenameSave = ['DepthProportion_' locID '_' num2str(DepthWindow(2)) 'm.mat'];
%         eval (['save -mat ' filenameSave ' PropDepths', ' PropMedianTemps', ' PropMeanTemps', ' PropTempSD', ' PropTempVar'])
    
    %% load the right one and plot
    
    filenameLoad = ['DepthProportion_CA_' num2str(DepthWindow(2)) 'm.mat']; % load CA
    eval (['load -mat ' filenameLoad])
    propD_CA = PropDepths;
    Tmedian_CA = PropMedianTemps;
    Tmean_CA = PropMeanTemps;
    Tsd_CA = PropTempSD;
    Tvar_CA = PropTempVar;
    
    filenameLoad = ['DepthProportion_GOC_' num2str(DepthWindow(2)) 'm.mat']; % load GOC
    eval (['load -mat ' filenameLoad])
    propD_GOC = PropDepths;
    Tmedian_GOC = PropMedianTemps;
    Tmean_GOC = PropMeanTemps;
    Tsd_GOC = PropTempSD;
    Tvar_GOC = PropTempVar;
    
    figure; hold on
    H(1) = plot(Tmean_CA, propD_CA, '.','markersize', 25);
    errorbar_x(Tmean_CA, propD_CA, Tsd_CA, '.')
    H(2) = plot(Tmean_GOC, propD_GOC, '.r', 'markersize', 25);
    errorbar_x(Tmean_GOC, propD_GOC, Tsd_GOC, '.')
    set(gca, 'YLim', [0 0.6], 'XLim', [10 30]);
    xlabel({['Mean temperature and SD at ' num2str(DepthWindow(1)) '-' num2str(DepthWindow(2)) 'm (°C)' ]},'FontSize', 14, 'FontWeight', 'bold')
    ylabel({['Proportion of 24-hour day at '  num2str(DepthWindow(1)) '-' num2str(DepthWindow(2)) 'm']},'FontSize', 14, 'FontWeight', 'bold')
    legend(H,['CCS, n=' num2str(length(Tmean_CA))],['GOC, n=' num2str(length(Tmean_GOC))],'Location','NorthEast')
    title({['Proportion of 24-hour day and temperatures at ' num2str(DepthWindow(1)) '-' num2str(DepthWindow(2)) 'm (' labDeploy ')']; ...
        '5 CCS tags, 6 GOC tags'}, 'FontSize', 14, 'FontWeight', 'bold');
    FigName = ['DepthProp_GOC_CCS_' num2str(DepthWindow(2)) 'm' labDeploy '.pdf']; % when run with just forage tags, label it.
    annotate_JS(Mfilename, gcf, FigName)
%         orient landscape % save
%         print('-dpdf', [dirBoth FigName]);
    
end


%%

disp('Completed calcDepthtimes.m')
% ===== EOF [calcDepthtimes.m] ======
