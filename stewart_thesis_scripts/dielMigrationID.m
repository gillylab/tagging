% dielMigrationID.m
% based off of findDielMigration.m by Ashley Booth, but is more lax: just
% takes any cast on either side of sunrise or sunset. 
%  
% Called From: 
% processSeriesData.m % J. Stewart
  
% Description: 
%  
% Outside Functions Called: 
%  
% AUTHOR    : J. Stewart jules32@gmail.com 
% DATE      : 26-Apr-2011 10:47:07  
% Revision  : 1.00  
% DEVELOPED : 7.9.0.529 (R2009b) OSX 
% FILENAME  : dielMigrationID.m 

Mfilename = mfilename; 

%% OLD way: Find Dusk Ascents and Dawn Descents
% dates = unique(floor(Seriesdate));
% 
% DepthDA2 = nan(length(DEPTH),1);
% DepthDD2 = nan(length(DEPTH),1);
% 
% buffer1 = datenum(0000,00,00,01,30,00);%within 1hr 30min of sunset
% buffer2 = datenum(0000,00,00,02,00,00);%within 1hr 30min of sunset
% 
% for dy = 1:length(DeployDays)  
%   %Dusk Ascents
%     b4sunset = dy + (sset - buffer1); %1hr 15min before sunset
%     aftsunset = dy + (sset + buffer1); %1hr 15min after sunset
%     if Seriesdate(1)<b4sunset %if dusk for day is after start of deployment    
%         startindex = Seriesdate(abs(Seriesdate-b4sunset) == min(abs(Seriesdate-b4sunset))); %find where buffer starts in day n
%         stopindex = Seriesdate(abs(Seriesdate-aftsunset) == min(abs(Seriesdate-aftsunset)));   
%    
%       % before sunset in relation to entire DEPTH array:
%         junkdepth = nan(stopindex,1); %create a blank array
%         junkdepth(startindex:stopindex) = DEPTH(startindex:stopindex); %fill array with only DEPTH data within sun buffer
%     
%         startDA = find(junkdepth == min(junkdepth)); %start index of DA = deepest point
%         stopDA = find(junkdepth == max(junkdepth)); %start index of DA = shallowest point
% 
%         DepthDA2(startDA:stopDA) = DEPTH(startDA:stopDA); %fill empty array with day(dy) with depths within DA
%     end
%     
%   %Dawn Descents
%     b4sunrise = dy + (srise - buffer2); %1hr 15min before sunrise
%     aftsunrise = dy + (srise + buffer2); %1hr 15min after sunrise
%     
%     if Seriesdate(1)<b4sunrise %if dusk for day is after start of deployment
%         startindex = Seriesdate(abs(Seriesdate-b4sunrise) == min(abs(Seriesdate-b4sunrise))); %find where buffer starts in day n
%         stopindex = Seriesdate(abs(Seriesdate-aftsunrise) == min(abs(Seriesdate-aftsunrise)));   
%     
%       % before sunrise in relation to entire DEPTH array:  
%         junkdepth = nan(stopindex,1); %create a blank array
%         junkdepth(startindex:stopindex) = DEPTH(startindex:stopindex); %fill array with only DEPTH data within sun buffer
%     
%         startDD = find(junkdepth == max(junkdepth)); %start index of DD = shallowest point
%         stopDD = find(junkdepth == min(junkdepth)); %start index of DD = deepest point
%         DepthDD2(startDD:stopDD) = DEPTH(startDD:stopDD); %fill empty array with day(dy) with depths within DD
%     end    
% end
% 
% % Plot day and night boxes with Dusk Ascents and Dawn Descents
% if exist('DayNight_boxes.m','file') %if the function is present the figure will have night boxes
%     [xb,yb] = DayNight_boxes(Seriesdate,sr/24,ss/24,[min(DEPTH),0],0); %will plot boxes
% %     fighandles(end+1) = figure(3);
%     set(gcf,'Position',[12   213   842   688])
% %     xlim = get(gca,'XLim'); ylim = get(gca,'YLim');
%     clf('reset');
%     fill(xb,yb,[0.8 0.8 0.8], 'EdgeColor','none')
%     hold on
%     h(1) = plot(Seriesdate,DEPTH,'k');
% %     plot(Seriesdate,DepthDD2,'-b');
% %     plot(Seriesdate,DepthDA2,'-r');
% %     h(2) = plot(Seriesdate,DepthDD,'-c.');
% %     h(3) = plot(Seriesdate,DepthDA,'-m.');
%     hold off
%     xlabel('Date')
%     ylabel('DEPTH [m]')
% %     title([name ' Dawn and Dusk Migrations'])
%     legend(h,'DEPTH','Dawn Descents','Dusk Ascents','Location','Southwest')
% %     set(gca,'XLim',xlim,'YLim',ylim)
%     datetick('x',6,'keeplimits')
% %     Annotate_plot([mfilename '.m'],gcf) %function bay A. Booth Apr 2009
% end


%% Find Dusk Ascents and Dawn Descents
% change in DEPTH over 30min is >100m w/o a chnge in DEPTH of 30m
% NOTE: uses negative depths
wiggle = 30; %can't wiggle by more than n meters
bigblock = 40*60;
timeblock = 10*60; %sec--size of blocks to look at 
dates = unique(floor(Seriesdate));

DepthDA = nan(length(DEPTH),1);
DepthDD = nan(length(DEPTH),1);

buffer1 = datenum(0000,00,00,02,00,00);%within 1hr 30min of sunset/rise
buffer2 = datenum(0000,00,00,02,00,00);%within 2hr of sunset/rise

for dy = min(dates):1:max(dates) %dates(13)
%%Dusk Ascents ---need to start from end and go back?
    b4sunset = dy + (sset - buffer1); 
    aftsunset = dy + (sset + buffer2); 
    jkindex = [];
    startDA = [];
    stopDA = [];
%     DAmeans=[];
    if(Seriesdate(1)<b4sunset) %if dusk for day is after start of deployment    
        startindex = find(abs(Seriesdate-b4sunset) == min(abs(Seriesdate-b4sunset))); %find where buffer starts in day n
        stopindex = find(abs(Seriesdate-aftsunset) == min(abs(Seriesdate-aftsunset)));   
      %to find index of deepest and shallowest point within 1hr 15min
      % before sunset in relation to entire DEPTH array:
        junkdepth = nan(stopindex,1); %create a blank array
        junkdepth(startindex:stopindex) = DEPTH(startindex:stopindex); %fill array with only DEPTH data within sun buffer
        startindex = find(junkdepth == min(junkdepth));
      %track through all depths within start and stop index
        jk = startindex;
        while jk < stopindex%jk = startindex:120:stopindex%start at begining of depths after nans
            ed = jk+bigblock; %30 min
            if ed > length(junkdepth)%make sure does not exceed matrix limit
                ed = length(junkdepth);
            end
          %this is a change in DEPTH over change in Seriesdate constraint
%                 DAmeans(end+1,1) = mean(junkdepth(jk:ed)); %mean of depths within 40min block
%                 DAmeans(end,2) = junkdepth(jk)+20;% start DEPTH
%                 DAmeans(end,3) = diff([mean(junkdepth(jk:ed)),junkdepth(jk)+20]);
%                 DAmeans(end,4) = mean(junkdepth(jk:ed))>junkdepth(jk)+20;
%                 DAmeans(end,5) = junkdepth(ed);
            if mean(junkdepth(jk:ed))>junkdepth(jk)+35 %find where mean change in DEPTH is > n m from start DEPTH
                jkindex = [jkindex,jk]; %keep track of past points
                %now look at finer scale for wiggles more than n meters
                for hg = jk:1:ed %$now cycle within 50 min block
                    tb = hg+timeblock;%$
                    if tb > ed %make sure does not exceed matrix limit
                        tb = ed;
                    end
                    block = ((junkdepth(hg:tb)-junkdepth(hg))<-wiggle);% take away current point from next n points and see if exceeds wiggle room
                    if sum(block)>0 %if there is a wiggle, find shallowest point before decrease by useing second derivative
                        sumblock = sum(block);
                        troughindex = find(junkdepth(1:tb)==max(junkdepth(jkindex(1):tb))); % find(diff(junkdepth(jk:tb))==0);
                        startDA = jkindex(1);%find(junkdepth(jkindex(1):jkindex(1)+bigblock) == min(junkdepth(jkindex(1):jkindex(1)+bigblock))); %jkindex(1) start index of DA = deepest point
                        stopDA = troughindex; %start index of DA = shallowest point
                        DepthDA(startDA:stopDA) = DEPTH(startDA:stopDA); %fill empty array with day(dy) with depths within DA
                        break
                    elseif (jk == stopindex)
                        startDA = jkindex(1); %start index of DA = deepest point
                        stopDA = hg; %start index of DA = shallowest point
                        DepthDA(startDA:stopDA) = DEPTH(startDA:stopDA); %fill empty array with day(dy) with depths within DA                        
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
                DepthDA(startDA:stopDA) = DEPTH(startDA:stopDA); %fill empty array with day(dy) with depths within DD
                break
            else
                jk =jk + 10*60; %create a larger search radius if not getting
            end
        end
    end
%%Dawn Descents ---need to start from end and go back?
    b4sunrise = dy + (srise - buffer2); %1hr 15min before sunrise
    aftsunrise = dy + (srise + buffer1); %1hr 15min after sunrise
    jkindex = [];
    startDD = [];
    stopDD = [];
    jk=[];
    if(Seriesdate(1)<b4sunrise) %if dawn for day is after start of deployment    
        startindex = find(abs(Seriesdate-b4sunrise) == min(abs(Seriesdate-b4sunrise))); %find where buffer starts in day n
        stopindex = find(abs(Seriesdate-aftsunrise) == min(abs(Seriesdate-aftsunrise)));   
      %to find index of deepest and shallowest point within 1hr 15min
      % before sunset in relation to entire DEPTH array:
        junkdepth = nan(stopindex,1); %create a blank array
        junkdepth(startindex:stopindex) = DEPTH(startindex:stopindex); %fill array with only DEPTH data within sun buffer
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
          %this is a change in DEPTH over change in Seriesdate constraint
            if mean(junkdepth(jk:ed))<junkdepth(jk)-30 %want the mean of the depths to be Xm deeper than the start DEPTH % sum((junkdepth(jk:ed)-junkdepth(jk))>100) %find where change in DEPTH is >100m in 30min
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
                        DepthDD(startDD:stopDD) = DEPTH(startDD:stopDD); %fill empty array with day(dy) with depths within DD
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
                DepthDD(startDD:stopDD) = DEPTH(startDD:stopDD); %fill empty array with day(dy) with depths within DD
                break
            else
                jk =jk+10*60;%create a larger search radius if not getting
            end
        end
    end
end
%load DivingAnalysis_DawnDusk 

% clear junkdepth wiggle timblock dates buffer1 buffer2 jkindex startDA stopDA startDD stopDD jk hg ed tb bigblock startindex stopindex b4sunrise aftsunrise b4sunset aftsunset

%% Plot day and night boxes with Dusk Ascents and Dawn Descents
if exist('DayNight_boxes.m','file') %if the function is present the figure will have night boxes
    [xb,yb] = DayNight_boxes(Seriesdate,sr/24,ss/24,[min(DEPTH),0],0); %will plot boxes
%     fighandles(end+1) = figure(3);
    set(gcf,'Position',[12   213   842   688])
%     xlim = get(gca,'XLim'); ylim = get(gca,'YLim');
    clf('reset');
    fill(xb,yb,[0.8 0.8 0.8], 'EdgeColor','none')
    hold on
    h(1) = plot(Seriesdate,DEPTH,'k');
%     plot(Seriesdate,DepthDD2,'-b');
%     plot(Seriesdate,DepthDA2,'-r');
    h(2) = plot(Seriesdate,DepthDD,'-c.');
    h(3) = plot(Seriesdate,DepthDA,'-m.');
    hold off
    xlabel('Date')
    ylabel('DEPTH [m]')
    title([tagNumS ' Dawn and Dusk Migrations'])
    legend(h,'DEPTH','Dawn Descents','Dusk Ascents','Location','Southwest')
%     set(gca,'XLim',xlim,'YLim',ylim)
    datetick('x',6,'keeplimits')
%     Annotate_plot([mfilename '.m'],gcf) %function bay A. Booth Apr 2009
end


%% 
  
disp('Completed dielMigrationID.m') 
% ===== EOF [dielMigrationID.m] ======  
