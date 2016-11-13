% compareTagsVertVel.m
%
% Called From:

% Description:%% calculate and plot vertical migrations based on sunrise
% and sunset.

%
% Outside Functions Called: compareTags.m
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 11-Nov-2011 16:03:44
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : compareTagsVertVel.m

Mfilename = mfilename;
wplot = 0;

% WHEN THIS WORKS, delete from CompareTagsFigs

%% date.

DateTagsVe = datevec(Rawtrix(:,1));

%% setup sunrise sunset % code from processSeriesData.m % J. Stewart

if TagsInCA
    lat = 36;% just one setting for all.
    lon = -122;
%     buffernum = 30; %buffernum = 1; % for 1 hour. But 30 mins is more reliable. smaller swath but less garbage.
%     buffer = datenum(0000,00,00,00,buffernum,00); % buffer = datenum(0000,00,00,1,00,00);
    
else
    lat = 27.5;
    lon = -111.5;
%     buffernum = 1;
%     buffer = datenum(0000,00,00,buffernum,00,00); % buffer = datenum(0000,00,00,1,00,00);
end

buffernum = 30; %buffernum = 1; % for 1 hour. But 30 mins is more reliable. smaller swath but less garbage.
    buffer = datenum(0000,00,00,00,buffernum,00); % buffer = datenum(0000,00,00,1,00,00);
   

%% run through for ascent and descent

V_DescentTrix_S = nan(30,length(TagIDs)); % will populate with velocities
V_AscentTrix_S = nan(30,length(TagIDs));


for u = 1:length(TagIDs)
    ulog = TagIDs(u) == Rawtrix(:,6);
    DateTagsV = datenum(0,0,0,DateTagsVe(ulog,4),DateTagsVe(ulog,5),DateTagsVe(ulog,6));
    DateTagsN = datenum(DateTagsV);
    
    % Get calculated sunrise/sunset times:
    if exist('sunrise_.m', 'file') %checks to see if W. Broenkow's files are present
        
        yr = floor(mean(DateTagsVe(ulog,1)));
        mon = floor(mean(DateTagsVe(ulog,2)));
        da = floor(mean(DateTagsVe(ulog,3))); %not exactly logical but sufficient for to find sunrise/sunset
        
        [sr,azimuth] = sunrise_(yr,mon,da,lat,lon*-1,1); %function by W. Broenkow of MLML
        srisea = (datevec(sr/24));
        srise = datenum(0,0,0,srisea(:,4), srisea(:,5),srisea(:,6)); % this way is correct, 17-Mar-2011
        
        [ss,azimuth] = sunset_(yr,mon,da,lat,lon*-1,1); %function by W. Broenkow of MLML
        sseta = (datevec(ss/24));
        sset = datenum(0,0,0,sseta(:,4), sseta(:,5),sseta(:,6));
        
    else %approximate
        disp('for accurate sunrise/sunset times download air_sea toolbox from http://woodshole.er.usgs.gov/operations/sea-mat/index.html')
        disp('or get sunrise_.m function by W. Broenkow of Moss Landing Marine Labs')
        disp('approximate sunrise/sunset hours: mar set=1:37GMT rise= 13:39GMT; oct set = 24:50, rise = 13.35; Jul set =2:24, rise = 12:45')
        if mean(Time_at_DepthData.TatDmon)<=6 | mean(Time_at_DepthData.TatDmon)>=10 %Winter
            sset =1;
            srise = 12;
        else % summer
            sset = 2;
            srise = 13;
        end
    end
    
    
    for w = 1:2
        
        % identify whether descent or ascent based on sunrise/sunset
        if w == 1 % descent
            iVert =  (DateTagsV>=(srise-buffer)) & (DateTagsV<(srise+buffer));
            wlabel = 'Descent';
            wcolor = [0 0.5020 0];
        elseif w == 2 % ascent
            iVert = (DateTagsV>=0) & (DateTagsV<(sset+buffer));
            wlabel = 'Ascent';
            wcolor = 'b';
        end
        
        % create vectors of each variable for ascent/descent
        TimeVert = Rawtrix(ulog,1); TimeVert(~iVert) = NaN;
        DepthVert = Rawtrix(ulog,2); DepthVert(~iVert) = NaN;
        TempVert = Rawtrix(ulog,3); TempVert(~iVert) = NaN;
        DeployDayVert = Rawtrix(ulog,5); DeployDayVert(~iVert) = NaN;
        
        % isolate vectors
        DepthVerti = DepthVert(iVert);
        TempVerti = TempVert(iVert);
        TimeVerti = TimeVert(iVert);
        DeployDayVerti = DeployDayVert(iVert);
        DeployDayUni = unique(DeployDayVerti);
        
        VertVelo_H = [];
        VertVelo_S = [];
        Dminmax = []; % for plotting/troubleshooting
        Tminmax = [];
        for r = 1:length(DeployDayUni)
            rlog = DeployDayVerti == DeployDayUni(r); % already takes ulog into account (DeployDay)
            
            Dtemp = DepthVerti(rlog); %depth
            Ttemp = TimeVerti(rlog); %time
            tM = Ttemp(find(Dtemp == max(Dtemp))); % in case multiple times at same depth
            tm = Ttemp(find(Dtemp == min(Dtemp)));
            
            if ~isempty(tm)
                Dminmax = [Dminmax; min(Dtemp) max(Dtemp)]; % identify delta depth
                Tminmax = [Tminmax; tm(1) tM(1)]; % identify delta time based on depths
            else
                Dminmax = [Dminmax; NaN NaN];
                Tminmax = [Tminmax; NaN NaN];
            end
        end
        if w == 2 % account for the problem that ascents cross DeployDay designations (sunset)
            Dminmax(:,1) = [Dminmax(2:end,1); NaN];
            Tminmax(:,1) = [Tminmax(2:end,1); NaN];
            Dminmax(:,2) = [Dminmax(1:end-2,2); Dminmax(end,2); NaN];
            Tminmax(:,2) = [Tminmax(1:end-2,2); Tminmax(end,2); NaN];
%             Dminmax(end,:) = []; % because that is a value for tomorrow which won't come
%             Tminmax(end,:) = [];
            %             VertVelo_H = (Dminmax(:,2)-Dminmax(:,1))./(Tminmax(:,2)-Tminmax(:,1));
        end
        
        a = abs(Tminmax(:,2)-Tminmax(:,1))*24;
        alog = a > 1.6; % incase the ascent didn't line up, will discard them.
        a(alog) = [];
        Dminmax(alog,:) = [];
        Tminmax(alog,:) = [];
            
        VertVelo_H = (Dminmax(:,2)-Dminmax(:,1))./a;
        VertVelo_S = VertVelo_H/(60*60);
        
        % populate matrices
        if w == 1
            V_DescentTrix_S(1:length(VertVelo_S),u) = VertVelo_S;
        elseif w == 2
            V_AscentTrix_S(1:length(VertVelo_S),u) = VertVelo_S;
        end
        
        % plot
        if wplot
            figure; hold on;
            set(gcf,'Position',[12   213   842   688])
            plot(Rawtrix(ulog,1)-datenum(0,0,0,8,0,0), Rawtrix(ulog,2), '.', 'color', [.8 .8 .8]);
            plot(TimeVerti-datenum(0,0,0,8,0,0),DepthVerti, '.', 'Color', wcolor, 'markersize', 8);
            plot(Tminmax(:,1)-datenum(0,0,0,8,0,0),Dminmax(:,1),'ro', 'markersize', 10)
            plot(Tminmax(:,2)-datenum(0,0,0,8,0,0),Dminmax(:,2),'ko', 'markersize', 10)
%             plot(Tminmax(end-1,1)-datenum(0,0,0,8,0,0),Dminmax(end-1,1),'ro', 'markersize', 20)
%                plot(Tminmax(end-1,2)-datenum(0,0,0,8,0,0),Dminmax(end-1,2),'ko', 'markersize', 20)
            %  
%             plot(Tminmax(end-1,1)-datenum(0,0,0,8,0,0),Dminmax(end-1,1),'ro', 'markersize', 20)
%                plot(Tminmax(end,2)-datenum(0,0,0,8,0,0),Dminmax(end,2),'ko', 'markersize', 20)
            %                         plot(([TimeVerti(a); TimeVerti(end)])-datenum(0,0,0,8,0,0), [DepthVerti(a); DepthVerti(end)], 'ro', 'markersize', 20)
            %                         plot(TimeVerti(a2+1)-datenum(0,0,0,8,0,0), DepthVerti(a2+1), 'ko', 'markersize', 20)
            ylabel('Depth [m]', 'FontSize',20, 'FontWeight', 'bold')
            xlabel('Time, PST', 'FontSize',20, 'FontWeight', 'bold')
            title([ locID '-' num2str(TagIDs(u)) ' ' wlabel ', ' num2str(buffernum) ' min buffer  '], 'FontSize',20, 'FontWeight', 'bold')
            datetick('x', 13);
            set(gca, 'YDir', 'reverse');
            grid on
            FigName = [locID '-' num2str(TagIDs(u)) '_Profile_' wlabel '.pdf'];
            annotate_JS(Mfilename, gcf, FigName)
            %                 orient landscape % save
            %                 print('-dpdf', [dirCompare FigName]);
        end
        
    end
end

%% save

% V_DescentTrix_S(V_DescentTrix_S < 0) = NaN; % get rid of ones that aren't ascents or descents-going the wrong way
% % V_DescentTrix_S(V_DescentTrix_S > 0) = NaN; % GOC--probs
%
if w == 1
    mean(nanmean(V_DescentTrix_S))
    std(nanstd(V_DescentTrix_S))
    median(nanmedian(V_DescentTrix_S))
elseif w == 2
    mean(nanmean(V_AscentTrix_S))
    std(nanstd(V_AscentTrix_S))
    median(nanmedian(V_AscentTrix_S))
end

mean(nanmean([V_DescentTrix_S; V_AscentTrix_S]))
std(nanstd([V_DescentTrix_S; V_AscentTrix_S]))
median(nanmedian([V_DescentTrix_S; V_AscentTrix_S]))

%%

disp('Completed compareTagsVertVel.m')
% ===== EOF [compareTagsVertVel.m] ======
