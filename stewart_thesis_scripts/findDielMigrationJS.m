function [DepthDA,DepthDD] = findDielMigrationJS(DepthTag,Seriesdate,Dmax,Dmin,ss,sr)
% findDielMigrationJS.m finds where squid does dawn decents and dusk ascents
% based on a DepthTag trace and a Seriesdate vector. This function is called by
% plotVertMigration. 
%  
% Sample Call: 
%  [DepthDA,DepthDD,fighandles] = findDielMigrationJS(DepthTag,Seriesdate,Dmax,Dmin,ss,sr,fighandles,name)
% Inputs: 
%  variables DepthTag,Seriesdate,Dmax,Dmin,sset,srise,fighandles from DivingAnalysis2.m
% Outputs: 
%  [DepthDA,DepthDD,fighandles
% Outside Functions Called: 
%  
% AUTHOR    : A. Booth ashley@boothfamily.com 
% DATE      : 19-Dec-2008 10:46:05, modified by J. Stewart 2 June 2011  
% Revision  : 1.00  
% DEVELOPED : 7.4.0.287 (R2007a) Windows XP 
% FILENAME  : findDielMigrationJS.m 

Mfilename = mfilename;
name = 'julie'

%%
srise = (sr/24)-floor(sr/24); %Matlab format
% sriseS = datestr(srise,13); %str format with only decimal(Seriesdate) data
% 
% [ss,azimuth] = sunset_(yr,mon,da,lat,lon*-1,1); %function by W. Broenkow of MLML
sset = (ss/24)-floor(ss/24);
% ssetS = datestr(sset,13);
% 

%% Find Dusk Ascents and Dawn Descents
% change in DepthTag over 30min is >100m w/o a chnge in DepthTag of 30m
% NOTE: uses negative depths
wiggle = 30; %can't wiggle by more than wiggle meters
bigblock = 40*60;
timeblock = 10*60; %sec--size of blocks to look at 
dates = unique(floor(Seriesdate));

DepthDA = nan(length(DepthTag),1);
DepthDD = nan(length(DepthTag),1);

buffer1 = datenum(0000,00,00,02,00,00);%within 1hr 30min of sunset/rise
buffer2 = datenum(0000,00,00,02,00,00);%within 2hr of sunset/rise

for dy = min(dates):1:max(dates) %dates(13)
%%Dusk Ascents ---need to start from end and go back?
    b4sunset = dy + (sset - buffer1); %1hr 15min before sunset
    aftsunset = dy + (sset + buffer2); %1hr 15min after sunset
    jkindex = [];
    startDA = [];
    stopDA = [];
%     DAmeans=[];
    if(Seriesdate(1)<b4sunset) %if dusk for day is after start of deployment    
        startindex = find(abs(Seriesdate-b4sunset) == min(abs(Seriesdate-b4sunset))); %find where buffer starts in day n
        stopindex = find(abs(Seriesdate-aftsunset) == min(abs(Seriesdate-aftsunset)));   
      %to find index of deepest and shallowest point within buffer
      % before sunset in relation to entire DepthTag array:
        junkdepth = nan(stopindex,1); %create a blank array
        junkdepth(startindex:stopindex) = DepthTag(startindex:stopindex); %fill array with only DepthTag data within sun buffer
        startindex = find(junkdepth == min(junkdepth));
      %track through all depths within start and stop index
        jk = startindex;
        while jk < stopindex%jk = startindex:120:stopindex%start at begining of depths after nans
            ed = jk+bigblock; %30 min
            if ed > length(junkdepth)%make sure does not exceed matrix limit
                ed = length(junkdepth);
            end
          %this is a change in DepthTag over change in Seriesdate constraint
%                 DAmeans(end+1,1) = mean(junkdepth(jk:ed)); %mean of depths within 40min block
%                 DAmeans(end,2) = junkdepth(jk)+20;% start DepthTag
%                 DAmeans(end,3) = diff([mean(junkdepth(jk:ed)),junkdepth(jk)+20]);
%                 DAmeans(end,4) = mean(junkdepth(jk:ed))>junkdepth(jk)+20;
%                 DAmeans(end,5) = junkdepth(ed);
            if mean(junkdepth(jk:ed))>junkdepth(jk)+35 %find where mean change in DepthTag is > n m from start DepthTag
                jkindex = [jkindex,jk]; %keep track of past points
                %now look at finer scale for wiggles more than n meters
                for hg = jk:1:ed %$now cycle within 50 min block
                    tb = hg+timeblock;%$
                    if tb > ed %make sure does not exceed matrix limit
                        tb = ed;
                    end
                    block = ((junkdepth(hg:tb)-junkdepth(hg))<-wiggle);% take away current point from next n points and see if exceeds wiggle room
                    if sum(block)>0 %if there is a wiggle, find shallowest point before decrease by using second derivative
                        sumblock = sum(block);
                        troughindex = find(junkdepth(1:tb)==max(junkdepth(jkindex(1):tb))); % find(diff(junkdepth(jk:tb))==0);
                        startDA = jkindex(1);%find(junkdepth(jkindex(1):jkindex(1)+bigblock) == min(junkdepth(jkindex(1):jkindex(1)+bigblock))); %jkindex(1) start index of DA = deepest point
                        stopDA = troughindex; %start index of DA = shallowest point
                        DepthDA(startDA:stopDA) = DepthTag(startDA:stopDA); %fill empty array with day(dy) with depths within DA
                        break
                    elseif (jk == stopindex)
                        startDA = jkindex(1); %start index of DA = deepest point
                        stopDA = hg; %start index of DA = shallowest point
                        DepthDA(startDA:stopDA) = DepthTag(startDA:stopDA); %fill empty array with day(dy) with depths within DA                        
                    end
                end
                if ~isempty(stopDA)
                    break
                end
                jk =jk + 2*60; %progress bigblock forward by n mins
            elseif ~isempty(jkindex) %if the squid does not acend as much as expected then get the index of the max within range of bigblock
                troughindex = find(junkdepth(1:ed)==max(junkdepth(jkindex(1):ed)));
                startDA = jkindex(1); %start index of DA = deepest point
                stopDA = troughindex; %start index of DA = shallowest point
                DepthDA(startDA:stopDA) = DepthTag(startDA:stopDA); %fill empty array with day(dy) with depths within DD
                break
            else
                jk =jk + 10*60; %create a larger search radius if not getting
            end
        end
    end
%%Dawn Descents ---need to start from end and go back?
    b4sunrise = dy + (srise - buffer2); 
    aftsunrise = dy + (srise + buffer1); 
    jkindex = [];
    startDD = [];
    stopDD = [];
    jk=[];
    if(Seriesdate(1)<b4sunrise) %if dawn for day is after start of deployment    
        startindex = find(abs(Seriesdate-b4sunrise) == min(abs(Seriesdate-b4sunrise))); %find where buffer starts in day n
        stopindex = find(abs(Seriesdate-aftsunrise) == min(abs(Seriesdate-aftsunrise)));   
      %to find index of deepest and shallowest point within 1hr 15min
      % before sunset in relation to entire DepthTag array:
        junkdepth = nan(stopindex,1); %create a blank array
        junkdepth(startindex:stopindex) = DepthTag(startindex:stopindex); %fill array with only DepthTag data within sun buffer
        startindex = find(junkdepth == max(junkdepth));% get shllowest point
     %% track through all depths within start and stop index
        jk = startindex;
        while jk < stopindex%for jk = startindex:120:stopindex%start at begining of depths after nans
%             if sum(jk == [startindex:540:stopindex])
%                 disp(num2str(jk))
%             end
            ed = jk + bigblock; %sec
            if ed > length(junkdepth)%make sure does not exceed matrix limit
                ed = length(junkdepth);
            end
          %this is a change in DepthTag over change in Seriesdate constraint
            if mean(junkdepth(jk:ed))<junkdepth(jk)-30 %want the mean of the depths to be Xm deeper than the start DepthTag % sum((junkdepth(jk:ed)-junkdepth(jk))>100) %find where change in DepthTag is >100m in 30min
                jkindex = [jkindex,jk]; %keep track of past points
                %now look at finer scale for wiggles more than n meters
                for hg = jk:1:ed %now cycle within big block to look for wiggles
                    tb = hg+timeblock;
                    if tb > ed %make sure does not exceed matrix limit
                        tb = ed;
                    end
                    block = ((junkdepth(hg:tb)-junkdepth(hg))>wiggle);% take away current point from next n points and see if exceeds wiggle room
                    if sum(block)>0 %if there is a wiggle, find deepest point before increase by useing second derivative
                        sumblock = sum(block);
                        peakindex = find(junkdepth(1:tb)==min(junkdepth(jkindex(1):tb))); % find the deepest point%find(diff(junkdepth(jk:tb))==0);
                        startDD = jkindex(1); %start index of DD = deepest point
                        stopDD = peakindex; %start index of DD = shallowest point
                        DepthDD(startDD:stopDD) = DepthTag(startDD:stopDD); %fill empty array with day(dy) with depths within DD
                        break
                    end
                end
                if ~isempty(startDD)%get out of loop if startDA has been filled
                    break
                end
                jk =jk+2*60; %progess search forward by n min
            elseif ~isempty(jkindex)%if the squid does not dive as much as expected then get the index of the min within range of bigblock
                peakindex = find(junkdepth(1:ed)==min(junkdepth(jkindex(1):ed)));
                startDD = jkindex(1); %start index of DD = shallowest point
                stopDD = peakindex; %start index of DD = deepest point
                DepthDD(startDD:stopDD) = DepthTag(startDD:stopDD); %fill empty array with day(dy) with depths within DD
                break
            else
                jk =jk+10*60;%create a larger search radius if not getting
            end
        end
    end
end

clear junkdepth wiggle timblock dates buffer1 buffer2 jkindex startDA stopDA startDD stopDD jk hg ed tb bigblock startindex stopindex b4sunrise aftsunrise b4sunset aftsunset

%% Plot day and night boxes with Dusk Ascents and Dawn Descents
if exist('DayNight_boxes.m','file') %if the function is present the figure will have night boxes
    [xb,yb] = DayNight_boxes(Seriesdate,sr/24,ss/24,[min(DepthTag),0],0); %will plot boxes
%     fighandles(end+1) = figure(3);
    set(gcf,'Position',[12   213   842   688])
%     xlim = get(gca,'XLim'); ylim = get(gca,'YLim');
    clf('reset');
    fill(xb,yb,[0.8 0.8 0.8], 'EdgeColor','none')
    hold on
    h(1) = plot(Seriesdate,DepthTag,'k');
%     plot(Seriesdate,DepthDD2,'-b');
%     plot(Seriesdate,DepthDA2,'-r');
    h(2) = plot(Seriesdate,DepthDD,'-c.');
    h(3) = plot(Seriesdate,DepthDA,'-m.');
%     hold off
    xlabel('Date')
    ylabel('DepthTag [m]')
    title([name ' Dawn and Dusk Migrations'])
    legend(h,'DepthTag','Dawn Descents','Dusk Ascents','Location','Southwest')
    set(gca, 'YDir', 'reverse');
%     set(gca,'XLim',xlim,'YLim',ylim)
    datetick('x',6,'keeplimits')
    Annotate_plot([mfilename '.m'],gcf) %function bay A. Booth Apr 2009
end

%% OLD way: Find Dusk Ascents and Dawn Descents. See A. Booth's script findDielMigration.m

%%

disp('Completed findDielMigrationJS.m')

% ===== EOF [findDielMigrationJS.m] ======
