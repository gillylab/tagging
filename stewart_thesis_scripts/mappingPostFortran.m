% mappingPostFortran.m
%
% Called From: habitatID_eachTag.m

% Description: this replaces radiusApproach_Tags.m: that is no longer
% necessary since the radiusApproach_RicardoPrep prepares data nicely for
% Fortran. This then is similar to R code 'RcodeRicardoSquidSBjunefortran.r': reads in the fortran output,
% maps it.
%
% Outside Functions Called:
% mappingPostFortran_dist.m % J. Stewart
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 26-Jun-2011 16:55:00
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : mappingPostFortran.m

Mfilename = mfilename;
% cd(dirFortran);

%% read in CalCOFI

% if TagHinke == 83051
%     manageTagData_CTD
%         cd(dirf); % important to have this here so that 83051 above won't change the directory!!!
% %     cd(dirFortran);
% end

%% plot maps

MarkerSizeID = [22 20 18 12];
sumdlog = [];
lsvectorinfo = [];
maxprobProbTrix = [];
for d = 1:length(probDayUni)
    figure; hold on
    plot(mapLon, mapLat, '.','color', [0.8 0.878 0.968], 'markersize', 30)
    dlog = probDayUni(d) == probDay; % called 'good' in R code
    plog = probProb(dlog) > .5; % called 'good' in R code
    sumdlog = [sumdlog sum(dlog)];
    probtemp = probProb(dlog);
    maxprobProb = max(probtemp); % minprobProbTrix. to find min, sort probtemp and look at next smallest that isn't 0
    maxprobProbTrix = [maxprobProbTrix; maxprobProb];
    abstrix = [];
    cmask = 1-probProb(dlog)/max(probProb(dlog));
    scatter(probLon(dlog),probLat(dlog), 100, cmask, 'filled','s') % now figure out how to display better
    colormap(gray(length(unique((cmask)))));
    text((min(mapLon)+range(mapLon)*.7), max(mapLat)-range(mapLat)*.05, ['black = ' num2str(maxprobProb,'%2.4f') ' '], 'color', 'k')
    
    % Ascents and Descents together
%     if TagHinke == 83051 % just do CalCOFI compare for tag 83051 tried with d not d-1
%         %         if d < length(probDayUni)
%         %             for g = 1:4
%         %                 plot(newLonA(lsindA(g,d)),newLatA(lsindA(g,d)), 'o', 'color', ccA(g,:), 'MarkerSize', 11, 'linewidth', 3)
%         %                 text((min(mapLon)+range(mapLon)*.7), max(mapLat)-g*range(mapLat)*.05,...
%         %                     [datestr(lsdateA(g,d), 'mmm-dd') ', rms=' num2str(lsvalA(g,d), '%2.2f') ' '], 'color', ccA(g,:))
%         %                 plot(newLonD(lsindD(g,d)),newLatD(lsindD(g,d)), 's', 'color', ccD(g,:), 'MarkerSize', 8, 'linewidth', 3)
%         %                 text((min(mapLon)+range(mapLon)*.7), (max(mapLat)-range(mapLat)*.2)-g*range(mapLon)*.05,...
%         %                     [datestr(lsdateD(g,d), 'mmm-dd') ', rms=' num2str(lsvalD(g,d), '%2.2f') ' '], 'color', ccD(g,:))
%         %             end
%         %         end
%         %         if d < length(probDayUni)
%         %             title({[num2str(TagHinke) ': Day ' num2str(d-1) ': ' datestr(HinkeDates(d),'mmm-dd') '  ']; ['Ascents: orange, descents: pink  ']}, 'FontSize',20, 'FontWeight', 'bold')
%         %         else
%         %             title([num2str(TagHinke) ': Day ' num2str(d-1) ': ' datestr(HinkeDates(end)+1,'mmm-dd') '  '], 'FontSize',20, 'FontWeight', 'bold')
%         %         end
%         if d == 6 || d == 7 || d == 16 || d == 17 % all 4 ascents. % don't do before day 5.
%             for g = 1:4
%                 aA = abs(HinkeDates(d)-lsdateA(g,d));
%                 plot(newLonA(lsindA(g,d)),newLatA(lsindA(g,d)), 'o', 'color', ccA(g,:), 'MarkerSize', 11, 'linewidth', 3)
%                 %                 text((min(mapLon)+range(mapLon)*.7), max(mapLat)-g*range(mapLat)*.05,...
%                 %                     [datestr(lsdateA(g,d), 'mmm-dd') ', rms=' num2str(lsvalA(g,d), '%2.2f') ' '], 'color', ccA(g,:))
%                 abstrix = [abstrix; aA];
%             end
%             title([num2str(TagHinke) ': Day ' num2str(d-1) ': ' datestr(HinkeDates(d),'mmm-dd') '  '], 'FontSize',20, 'FontWeight', 'bold')
%             sourceid = [0; 0; 0; 0]; % 0= ascent
%             lsvectorinfo = [lsvectorinfo; lsvalA(1:4,d) sourceid abstrix];
%             
%         elseif  d == 9 || d == 10 || d == 12 % all 4 descents . don't start till day 5. d == 2 || d == 4 || d == 5 ||
%             for g = 1:4
%                 aD = abs(HinkeDates(d)-lsdateD(g,d));
%                 plot(newLonD(lsindD(g,d)),newLatD(lsindD(g,d)), 's', 'color', ccD(g,:), 'MarkerSize', 8, 'linewidth', 3)
%                 %                 text((min(mapLon)+range(mapLon)*.7), (max(mapLat)-range(mapLat)*.2)-g*range(mapLon)*.05,...
%                 %                     [datestr(lsdateD(g,d), 'mmm-dd') ', rms=' num2str(lsvalD(g,d), '%2.2f') ' '], 'color', ccD(g,:)) %delete text when know it's right
%                 abstrix = [abstrix; aD];
%             end
%             title([num2str(TagHinke) ': Day ' num2str(d-1) ': ' datestr(HinkeDates(d),'mmm-dd') '  '], 'FontSize',20, 'FontWeight', 'bold')
%             sourceid = [1; 1; 1; 1]; % 1= descent
%             lsvectorinfo = [lsvectorinfo; lsvalD(1:4,d) sourceid abstrix];
%             
%         elseif d == 11 || d == 13 % a mix: 3 ascents and 1 descent % THIS IS DAY 10!!
%             for g = 1:3
%                 aA = abs(HinkeDates(d)-lsdateA(g,d));
%                 plot(newLonA(lsindA(g,d)),newLatA(lsindA(g,d)), 'o', 'color', ccA(g,:), 'MarkerSize', 11, 'linewidth', 3)
%                 abstrix = [abstrix; aA];       
%             end
%             plot(newLonD(lsindD(1,d)),newLatD(lsindD(1,d)), 's', 'color', ccD(g,:), 'MarkerSize', 8, 'linewidth', 3)
%             aD = abs(HinkeDates(d)-lsdateD(1,d));
%             abstrix = [abstrix; aD];
%             title([num2str(TagHinke) ': Day ' num2str(d-1) ': ' datestr(HinkeDates(d),'mmm-dd') '  '], 'FontSize',20, 'FontWeight', 'bold')
%             sourceid = [0; 0; 0; 1]; % 0= ascent
%             lsvectorinfo = [lsvectorinfo; [lsvalA(1:3,d); lsvalD(1,d)] sourceid abstrix];
%         elseif d == 8 ||d == 14 % a mix: 3 descents and 1 ascent
%             for g = 1:3
%                 aD = abs(HinkeDates(d)-lsdateD(g,d));
%                 plot(newLonD(lsindD(g,d)),newLatD(lsindD(g,d)), 'o', 'color', ccD(g,:), 'MarkerSize', 11, 'linewidth', 3)
%                 abstrix = [abstrix; aD];
%             end
%             plot(newLonA(lsindA(1,d)),newLatA(lsindA(1,d)), 's', 'color', ccA(g,:), 'MarkerSize', 8, 'linewidth', 3)
%             aA = abs(HinkeDates(d)-lsdateA(1,d));
%             abstrix = [abstrix; aA];
%             title([num2str(TagHinke) ': Day ' num2str(d-1) ': ' datestr(HinkeDates(d),'mmm-dd') '  '], 'FontSize',20, 'FontWeight', 'bold')
%             sourceid = [1; 1; 1; 0]; % 0= ascent
%             lsvectorinfo = [lsvectorinfo; [lsvalD(1:3,d); lsvalA(1,d)] sourceid abstrix];
%             
%         elseif d == 15 % a mix: 2 and 2.
%             for g = 1:2
%                 aA = abs(HinkeDates(d)-lsdateA(g,d));
%                 aD = abs(HinkeDates(d)-lsdateD(g,d));
%                 plot(newLonA(lsindA(g,d)),newLatA(lsindA(g,d)), 's', 'color', ccA(g,:), 'MarkerSize', 8, 'linewidth', 3)
%                 plot(newLonD(lsindD(g,d)),newLatD(lsindD(g,d)), 'o', 'color', ccD(g,:), 'MarkerSize', 11, 'linewidth', 3)
%                 abstrix = [abstrix; aA; aD];
%             end
%             title([num2str(TagHinke) ': Day ' num2str(d-1) ': ' datestr(HinkeDates(d),'mmm-dd') '  '], 'FontSize',20, 'FontWeight', 'bold')
%             sourceid = [0; 1; 0; 1]; % 0= ascent
%             lsvectorinfo = [lsvectorinfo; [lsvalD(1:3,d); lsvalA(1,d)] sourceid abstrix];
%         else
%             title([num2str(TagHinke) ': Day ' num2str(d-1) ': ' datestr(HinkeDates(d),'mmm-dd') '  '], 'FontSize',20, 'FontWeight', 'bold')
%         end
%         
%     else
%         title([num2str(TagHinke) ': Day ' num2str(d-1) ': ' datestr(HinkeDates(d),'mmm-dd') '  '], 'FontSize',20, 'FontWeight', 'bold')
%     end
    
    plot(LocTagCA_Hinke(:,1), LocTagCA_Hinke(:,2),'x', 'markersize', 10, 'linewidth', 3)
    FigName = [num2str(TagHinke) '_mapPostFortranAD_Day' num2str(d-1) '.pdf'];
    annotate_JS(Mfilename, gcf, FigName);
    orient landscape % save
%     print('-dpdf', [dirFortran FigName]);
end

% if TagHinke == 83051
%     lsvectorinfo % for arcgis
%     maxprobProbTrix
%     
%     % for day 10
%     aA = abs(HinkeDates(d)-lsdateA(g,d));
%     aD = abs(HinkeDates(d)-lsdateD(g,d));
%     plot(newLonA(lsindA(g,d)),newLatA(lsindA(g,d)), 's', 'color', ccA(g,:), 'MarkerSize', 8, 'linewidth', 3)
%     plot(newLonD(lsindD(g,d)),newLatD(lsindD(g,d)), 'o', 'color', ccD(g,:), 'MarkerSize', 11, 'linewidth', 3)
%     abstrix = [abstrix; aA; aD];
% end
%% % distance between each best/worst points between days--not super interesting. see mappingPostFortran_dist.m instead. 

DistTrix = nan(4, length(probDayUni)-1); % set up for min, mean, max distance for each day
% DistTrix2 = nan(3, length(probDayUni)-1);
% n2=2;
counthighest = [];
ProbTrix = [probProb probLat probLon];
for d = 1:length(probDayUni)
    if d < length(probDayUni)
        % id today
        dlog = probDayUni(d) == probDay;
        probtemp = probProb(dlog);
        maxprobProb = max(probtemp);
        x = probLon(dlog);
        y = probLat(dlog);
        mlog = probtemp == maxprobProb;
        mdex = find(probtemp == maxprobProb);% can erase afters
        xm = x(mlog);
        ym = y(mlog);
        
        % id tomorrow
        dlog2 = probDayUni(d+1) == probDay;
        probtemp2 = probProb(dlog2);
        maxprobProb2 = max(probtemp2);
        x2 = probLon(dlog2);
        y2 = probLat(dlog2);
        mlog2 = probtemp2 == maxprobProb2;
        mdex2 = find(probtemp2 == maxprobProb2); % can erase afters
        xm2 = x2(mlog2);
        ym2 = y2(mlog2);
        counthighest = [counthighest [length(xm); length(xm2)]];
        
        % distance between best positions
        disttrix = [];
        for i = 1:length(xm)
            [dist_km] = distance_haversine(xm(i), ym(i), xm2, ym2);
            disttrix = [disttrix; dist_km];
        end
        if size(disttrix,2) == 1
            DistTrix(:,d) = [min(disttrix); mode(disttrix); mean(disttrix); max(disttrix)];
        else
            DistTrix(:,d) = [min(min(disttrix)); mode(mode(disttrix)); mean(mean(disttrix)); max(max(disttrix))];
        end
        
        % distance between worst positions: NOT INTERESTING
    end
end
% DistTrix_km = distdim(DistTrix,'deg','km', 'earth');

% plot
figure; hold on
plot(1:length(probDayUni)-1,DistTrix(4,:), '-', 'color', [0.8 0.8 0.8], 'linewidth', 1.5)
plot(1:length(probDayUni)-1,DistTrix(3,:), 'k--', 'linewidth', 1.5)
plot(1:length(probDayUni)-1,DistTrix(2,:), 'k', 'linewidth', 3)
plot(1:length(probDayUni)-1,DistTrix(1,:), '--', 'color', [0.8 0.8 0.8], 'linewidth', 1.5)
legend('max', 'mean','mode', 'min');
xlabel('Day ','FontSize', 14, 'FontWeight', 'bold')
ylabel('distance traveled (km)  ', 'FontSize', 14, 'FontWeight', 'bold')
title('Distances traveled each day (highest probability positions only)', 'FontSize', 20, 'FontWeight', 'bold')
FigName = [num2str(TagHinke) '_dist_traveled.pdf'];
annotate_JS(Mfilename, gcf, FigName);
% orient landscape % save
% print('-dpdf', [dirFortran FigName]);

%% plot timeseries figures of OLS v. time

MarkerSizeID = [12 10 8 2];
nn = 4;

countProb1 = []; % count for how many pixels with probability==1
for d = 1:length(probDayUni)
    dlog = probDayUni(d) == probDay; % called 'good' in R code
    plog = probProb(dlog) > 0; % called 'good' in R code
    countProb1 = [countProb1; sum(plog)]; % count
end

figure; hold on
for z = 1:length(probDayUni)
    zlog = probDayUni(z) == probDay; % called 'good' in R code
    probTemp = probProb(zlog);
    %         plot(repmat(z,length(probTemp),1),probTemp, '.');
    %         if z>1
    %             for g = 1:nn
    %
    %                 % do all this labeling with z-1 to be consistent with the tagging day = 0
    %                 aA = abs(HinkeDates(z)-lsdateA(g,z-1)); % ascents: orange
    %                 if round(aA) == 0
    %                     plot(z-1,lsvalA(g,z-1), '.', 'color', ccA(g,:), 'MarkerSize', 35)
    %                 elseif round(aA) <= nn && round(aA) ~= 0
    %                     plot(z-1,lsvalA(g,z-1), 'o', 'color', ccA(g,:), 'MarkerSize', MarkerSizeID(round(aA)), 'linewidth', 1.5)
    %                 else
    %                     plot(z-1,lsvalA(g,z-1), 'o', 'color', ccA(g,:), 'MarkerSize', MarkerSizeID(nn), 'linewidth', 1)
    %                 end
    %
    %                 aD = abs(HinkeDates(z)-lsdateD(g,z-1)); % descents: pink
    %                 if round(aD) == 0
    %                     plot(z-1,lsvalD(g,z-1), '.', 'color', ccD(g,:), 'MarkerSize', 35)
    %                 elseif round(aD) <= nn && round(aD) ~= 0
    %                     plot(z-1,lsvalD(g,z-1), 'o', 'color', ccD(g,:), 'MarkerSize', MarkerSizeID(round(aD)), 'linewidth', 1.5)
    %                 else
    %                     plot(z-1,lsvalA(g,z-1), 'o', 'color', ccA(g,:), 'MarkerSize', MarkerSizeID(nn), 'linewidth', 1)
    %                 end
    %             end
    %         end
    
    %         if TagHinke == 83051 % just do CalCOFI compare for tag 83051
    %         for g = 1:nn
    %
    %             % do all this labeling with z-1 to be consistent with the tagging day = 0
    %             aA = abs(HinkeDates(z)-lsdateA(g,z)); % ascents: orange
    %             if round(aA) == 0
    %                 plot(z-1,lsvalA(g,z), '.', 'color', ccA(g,:), 'MarkerSize', 35)
    %             elseif round(aA) <= nn && round(aA) ~= 0
    %                 plot(z-1,lsvalA(g,z), 'o', 'color', ccA(g,:), 'MarkerSize', MarkerSizeID(round(aA)), 'linewidth', 1.5)
    %             else
    %                 plot(z-1,lsvalA(g,z), 'o', 'color', ccA(g,:), 'MarkerSize', MarkerSizeID(nn), 'linewidth', 1)
    %             end
    %
    %             aD = abs(HinkeDates(z)-lsdateD(g,z)); % descents: pink
    %             if round(aD) == 0
    %                 plot(z-1,lsvalD(g,z), '.', 'color', ccD(g,:), 'MarkerSize', 35)
    %             elseif round(aD) <= nn && round(aD) ~= 0
    %                 plot(z-1,lsvalD(g,z), 'o', 'color', ccD(g,:), 'MarkerSize', MarkerSizeID(round(aD)), 'linewidth', 1.5)
    %             else
    %                 plot(z-1,lsvalA(g,z), 'o', 'color', ccA(g,:), 'MarkerSize', MarkerSizeID(nn), 'linewidth', 1)
    %             end
    %         end
    %         end
end
plot(0:length(probDayUni)-1, countProb1/sum(msk), '-', 'color', [0.8 0.8 .8], 'linewidth', 2)
xlabel('Day ','FontSize', 14, 'FontWeight', 'bold')
ylabel({'proportion of non-zero probability tagging area: gray'}, 'FontSize', 14, 'FontWeight', 'bold') % 'RMS of CalCOFI:tag 5 (orange: ascents, pink: descents) '; 
set(gca, 'fontsize', 12, 'fontweight', 'bold')
title({[num2str(TagHinke) ' Root mean squares (RMS) of CalCOFI v. tag-sampled temperature profiles ']; ...
    'and proportion of tagging area (ocean only) with non-zero probability values'}, 'fontsize', 16, 'fontweight', 'bold');
FigName = [num2str(TagHinke) '_RMS_timeseries.pdf'];
annotate_JS(Mfilename, gcf, FigName);
%     orient landscape % save
%     print('-dpdf', [dirFortran FigName]);

%% prepare for arcGIS ColorMaps

a = sort(probProb)
head(a(a > 0))


%%

disp('Completed mappingPostFortran.m')
% ===== EOF [mappingPostFortran.m] ======
