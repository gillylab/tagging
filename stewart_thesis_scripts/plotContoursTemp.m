% plotContoursTemp.m
%  based on contour_temp_profile_v2 by A. Carlisle. Called from
%  processSeriesData.m
%
%
% Outside Functions Called:
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 18-Mar-2011 16:31:37
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : plotContoursTemp.m

Mfilename = mfilename;

%%

%Variables you need to set to get particular results/output
temptype = 4; %4 = stalk TempTag (column 2), 3 = int TempTag (col 3), for use with archival record
DepthTagincr = 5; %size of DepthTag bins
gridinterp = 2; %set to 1 if use old griddata function, 2 if use newer TriScatteredInterp
addcontour = 0; %set to 1 if want TempTag contours on final figure, 0 if otherwise
median_mode_Z = 1; %use 1 if want to use median DepthTag, 0 if modal DepthTag
% depthLim = 300;  %set this to DepthTag you want the maximum DepthTag in the plot to be

depthLim=max(DepthTag);
tempLim=max(TempTag);
XI = fix(min(Seriesdate)):1:ceil(max(Seriesdate)); % 1 is 1 day. datenum(0,0,0,12,0,0)
YI = min(DepthTag):DepthTagincr:max(DepthTag);

a = isnan(Seriesdate);
b = isnan(DepthTag);
c = isnan(TempTag);
abc = a+b+c;
abclog = abc >= 1;


if gridinterp == 1  %old way, slower
    [xx,yy,zz] = griddata(Seriesdate,DepthTag,TempTag, XI, YI');  %xx-date, yy-DepthTag, zz-TempTag, forms grid with 3 layers
elseif gridinterp ==2  %updates function, does same as above but faster
    F = TriScatteredInterp(Seriesdate(~abclog),DepthTag(~abclog),TempTag(~abclog),'nearest');
%  F = TriScatteredInterp(Seriesdate(~abclog),DepthTag(~abclog));
    [xx,yy]=meshgrid(XI,YI);
    zz = F(xx,yy);
end

%find average DepthTag during day/night/overall day. Output is dayindex, a
%vector of 0's, 1's designating day or night

% [dayindex] = day_night_assignment_v2(light,DepthTag,Seriesdate);

%calculates daily mean/median DepthTag or TempTag, first parameter is data you
%want summarized, Seriesdate is the time vector, and need to enter a string
%with either 'mean' or 'median' depending on what stat you want. Daymeans
%is a 3xnumb days array, col 1 is mean Z at night, 2 is mean Z at day, 3 is
%overall mean.

% [xvec,yvec,daymeans,minZ,maxZ] = calcTSmeans(DepthTag,'DepthTag',Seriesdate,'mean',dayindex);

%get max DepthTag and set mask
date = (fix(Seriesdate));     %rounds all the dates down
dates = unique (date);      %picks out each unique date, so get subset of data w/ one row/day
% dates = XI;
offset= min(DepthTag)- 0;      %gets offset for DepthTag so that its min is 0 and no negaties

clear M;

daymax(1:length(dates))=nan;
daymax=daymax';
for y = 1:length(dates)
    daymax(y) = max(DepthTag(fix(Seriesdate) == dates(y)));
end

%% fill in single nans with mean of adjacent values first columns then rows
zzsub = zz;
for x = 1  : size(zz,2)
    ind = find(isnan(zz(:,x)));
    zzsub(ind(ind>1 & ind<size(zzsub,1)),x) = ...
        ( zzsub(ind(ind>1 & ind<size(zzsub,1))+1,x) + ...
        zzsub(ind(ind>1& ind<size(zzsub,1))-1,x))/2 ;
end

for x = 1  : size(zz,1)
    ind = find(isnan(zz(x,:)));
    zzsub(x,ind(ind>1 & ind<size(zzsub,2))) = ...
        ( zzsub(x,ind(ind>1 & ind<size(zzsub,2))+1) + ...
        zzsub(x,ind(ind>1& ind<size(zzsub,2))-1))/2 ;
end

zz = zzsub;
zzz=zz;

%pads 2 duplicates of first/last col and first/last row to compensate for
%smoothing artifacts
padnumb = 2;
sizemat=size(zzz);
zzz1(1:sizemat(1)+(2*padnumb),1:sizemat(2)+(2*padnumb))=NaN;
zzz1(padnumb+1:sizemat(1)+padnumb,padnumb+1:sizemat(2)+padnumb)=zzz;
for ab = padnumb:-1:1
    zzz1(ab,:)=zzz1(padnumb+1,:);
    zzz1(end-(ab-1),:)=zzz1(end-padnumb,:);
    zzz1(:,ab)=zzz1(:,padnumb+1);
    zzz1(:,end-(ab-1))=zzz1(:,end-padnumb);
end



%%smooth using conv2
F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
zzz1 = conv2(zzz1,F,'same');
zzz1 = conv2(zzz1,F,'same');
zzz1 = conv2(zzz1,F,'same');

%clip off duplicate first/last cols/rows
zzz1(1:padnumb,:)=[];
zzz1((end-(padnumb-1)):end,:)=[];
zzz1(:,1:padnumb)=[];
zzz1(:,(end-(padnumb-1)):end)=[];

zzz = zzz1;
clear zzz1;


%fill in NaN values at top of matrix so no blank areas at top of figure
zzzSize = size(zzz);
for ii = 1:zzzSize(2)
    znans = find(isnan(zzz(:,ii)));
    zdiff = diff(znans);
    zdiff2 =max(zdiff);
    if zdiff2>1
        zdiff2 = find(zdiff == zdiff2);
        zzz(1:zdiff2,ii) = zzz(zdiff2+1,ii);
    end
end


%% USES DAILY MAX TO CLIP
yy(:,:) = floor(yy(:,:));
%clip by max DepthTag of day
for aaa = 1:zzzSize(2)-1;  %-1 b/c padded extra day onto end for XI
    daymax(aaa) = ceil(daymax(aaa));
    tempz = daymax(aaa)/DepthTagincr;
    tempz = ceil(tempz)+1;
    daymaxZZZindex(aaa) = tempz;
end


for iii = 1:zzzSize(2)-1;  %-1 b/c padded extra day onto end for XI
    zzz(daymaxZZZindex(iii):end,iii)=NaN; %set all cells deeper than days max DepthTag as NaN
end

zeroIdx = find(yy(:,1)==0);
zzz(1:zeroIdx,:)=NaN; %set all values that are above 0 (above the surface) as NaN
% zzz(:,1) = zzz(:,2); %first column is all NaN's b/c of contouring, so set it to be equal to next days


%% plot

% just the interp
figure
scrsz = get(0,'ScreenSize');
% figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/4]);  %change size of figure to be
pcolor(xx,yy,zzz)
shading flat  %shading interp to make it smooth, but this can clip somewhat, so used flat to troubleshoot
% shading interp
set(gca,'YDir','reverse');
colorbar

datetick
ylim ([-1 depthLim])
xlim ([min(date) max(date)+1])  %use this to set xaxis to bound data set
%     xlim([184 486]);
hold on
if addcontour == 1
    [CS H] = contour(xx,yy,zzz,[2:2:30], 'k');
end
yvals = 0:100:depthLim;
set(gca,'YTick', yvals);
xlabel ('Date','fontsize', 16);
ylabel ('Depth (m)','fontsize', 16);
title([tagName ' Contour Depth Profile'], 'FontSize',20, 'FontWeight', 'bold')
FigName = [tagName '_Profile_ContourTemp'];
annotate_JS(Mfilename, gcf, FigName);
orient landscape % save
print('-dpdf', [FigName '.pdf']);
print('-dpng', [FigName '.png']);
print('-djpeg', [FigName '.jpg']);


% and with tag profile overlay
figure
scrsz = get(0,'ScreenSize');
pcolor(xx,yy,zzz)
shading flat  %shading interp to make it smooth, but this can clip somewhat, so used flat to troubleshoot
% shading interp
set(gca,'YDir','reverse');
colorbar
datetick
ylim ([-1 depthLim])
xlim ([min(date) max(date)+1])  %use this to set xaxis to bound data set
hold on
yvals = 0:100:depthLim;
set(gca,'YTick', yvals);
xlabel ('Date','fontsize', 16);
ylabel ('Depth (m)','fontsize', 16);
title({[tagName ' Contour Depth Profile']; ['Black = night; Gray = Day']}, 'FontSize',20, 'FontWeight', 'bold')
nightcol = [0.494117647409439 0.494117647409439 0.494117647409439];
hold on
plot(Seriesdate(iNight),DepthTag(iNight), '.k', 'MarkerSize', 15);
plot(Seriesdate(iDay),DepthTag(iDay), '.', 'color', nightcol, 'MarkerSize',15);
FigName = [tagName '_Profile_ContourTempTag'];
annotate_JS(Mfilename, gcf, FigName);
orient landscape % save
print('-dpdf', [FigName '.pdf']);
print('-dpng', [FigName '.png']);
print('-djpeg', [FigName '.jpg']);

%%

disp('Completed plotContoursTemp.m')
% ===== EOF [plotContoursTemp.m] ======
