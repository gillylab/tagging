% make_map_stewart.m

load mmusyl_marlin_all.mat % make this work for me!

% add my path of trusty ancillary matlab scripts
path(path,'~/mfiles');

% set upp map boundaries
xrange = max(lon) - min(lon);
yrange = max(lat) - min(lat);
xmin=min(lon)-.05*xrange;
xmax=max(lon)+.05*xrange;
ymin=min(lat)-.05*yrange;
ymax=max(lat)+.05*yrange;

% figure out IDs of unique beasties
idlist = unique(ptt);

% start plot
figure
m_proj('Mercator','longi',[xmin xmax],'lati',[ymin ymax]);
hold on
% add a rough coast
m_gshhs_c('patch',[.7 .7 .7]);

% This is where I could have it look up the Program Report files along
% with pop up lat/lon and have it plot it all
% add tracks one critter at a time
% for i = 1:length(idlist);
% 
%   % find all positions for the given beast
%   ind = find(ptt==idlist(i));
% 
%   % plot error bars based on the state-space 95% CI
%   xtmp = lon(ind);
%   ytmp = lat(ind);
%   ubx = hilon(ind);
%   lbx = lolon(ind);
%   uby = hilat(ind);
%   lby = lolat(ind);
% 
%   for j = 1:length(ind);
%    
%     % make horizontal line
%     h5 = m_plot([lbx(j) ubx(j)],[ytmp(j) ytmp(j)],'k');
%     h6 = m_plot([xtmp(j) xtmp(j)],[lby(j) uby(j)],'k');
%    
%     % make them light yet noticeable
%     set(h5,'Color',[.7 .7 .7]);
%     set(h6,'Color',[.7 .7 .7]);
%   end
% 
%   % plot trajectory as cyan line
%   m_plot(lon(ind),lat(ind),'c');
% 
% end



% add title
title('Marlins - postions with 95% CI used for data extractions')

% finalize plot
m_grid

% print it
print -dpng musyl_marlin_all_map.png
