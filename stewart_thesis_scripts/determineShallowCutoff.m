% determineShallowCutoff.m
%
% called from habitatID_eachTag.m
% this goes goes through the DailyShallow mat files and if there is enough
% data at that cutoff, it puts it into a matrix. Otherwise, it loads the
% next deepest for that day.
%
% Currently, minimum amount of data (nMin) = 5. Cutoffs go from 10-15-20
%
% Outside Functions Called:
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 22-Mar-2011 11:07:50
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : determineShallowCutoff.m

% *********TO BEGIN********
% find and replace all for the tag number you'd like. It will load and save everything properly
tagNum = 83052;

% Mfilename = mfilename; % no figs, no need to label
nMin = 5; % what is the fewest passable number of temp values

% to ensure all days are represented, first define uniDex with a deeper cutoff
load -mat DailyShallow_83052_50.mat % from tagDailyHistos_Shallow.m
uniDex = unique(ShallowDailyTemp(:,1));

%% first pass: the shallowest cutoff. 10m.

load -mat DailyShallow_83052_10.mat % from tagDailyHistos_Shallow.m
co = 10; %cutoff

uniDexCount = [];
ShallowDailyTempBest1 = [];
for i = 1:length(uniDex)
    ilog = uniDex(i) == ShallowDailyTemp(:,1);
    uniDexCount = [uniDexCount sum(ilog)]; % just to get a count
    if sum(ilog) >= nMin
        ShallowDailyTempBest1 = [ShallowDailyTempBest1; ShallowDailyTemp(ilog,:) repmat(co,sum(ilog),1)];
    end
end

uniDexBest1 = unique(ShallowDailyTempBest1(:,1));

%% second pass: 15m

if length(uniDexBest1) < length(uniDex)
    toCompleteDex2 = []; % identify which of uniDex are still to be completed
    for i = 1:length(uniDex)
        ilog = uniDexBest1 == uniDex(i);
        if sum(ilog) <= 0 % if this is not in uniDexBest1
            toCompleteDex2 = [toCompleteDex2; uniDex(i)];
        end
    end
    
    load -mat DailyShallow_83052_15.mat
    co = 15;
    
    uniDexCount2 = [];
    ShallowDailyTempBest2 = [];
    for i = 1:length(toCompleteDex2)
        ilog = toCompleteDex2(i) == ShallowDailyTemp(:,1);
        uniDexCount2 = [uniDexCount2 sum(ilog)];
        if sum(ilog) >= nMin
            ShallowDailyTempBest2 = [ShallowDailyTempBest2; ShallowDailyTemp(ilog,:) repmat(co,sum(ilog),1)];
        end
    end
    a = isempty(ShallowDailyTempBest2);
    if a == 1 % troubleshoot: if the 15m increment doesn't add anything
        ShallowDailyTempBest2 = ShallowDailyTempBest1;
        uniDexBest12 = [uniDexBest1];
    else
        uniDexBest2 = unique(ShallowDailyTempBest2(:,1));
        uniDexBest12 = [uniDexBest2; uniDexBest1];
    end
    if length(uniDexBest12) == length(uniDex)
        ShallowDailyTempBest = [ShallowDailyTempBest1; ShallowDailyTempBest2];
        save -mat DailyShallow_83052_Best.mat ShallowDailyTempBest DepthTagMax
    end
elseif length(uniDexBest1) == length(uniDex)
    ShallowDailyTempBest = ShallowDailyTempBest1;
    save -mat DailyShallow_83052_Best.mat ShallowDailyTempBest DepthTagMax
end


% combine the two

%% Third pass: 20m

if length(uniDexBest12) < length(uniDex)
    toCompleteDex3 = []; % figure out which of uniDex are still to be completed
    for i = 1:length(uniDex)
        ilog = uniDexBest12 == uniDex(i);
        if sum(ilog) <= 0 % if this is not in uniDexBest3
            toCompleteDex3 = [toCompleteDex3; uniDex(i)];
        end
    end
    
    load -mat DailyShallow_83052_20.mat
    co = 20;
    
    uniDexCount3 = [];
    ShallowDailyTempBest3 = [];
    for i = 1:length(toCompleteDex3)
        ilog = toCompleteDex3(i) == ShallowDailyTemp(:,1);
        uniDexCount3 = [uniDexCount3 sum(ilog)];
        if sum(ilog) >= nMin
            ShallowDailyTempBest3 = [ShallowDailyTempBest3; ShallowDailyTemp(ilog,:) repmat(co,sum(ilog),1)];
        end
    end
    a = isempty(ShallowDailyTempBest3);
    if a == 1 % troubleshoot: if the 15m increment doesn't add anything
        ShallowDailyTempBest3 = ShallowDailyTempBest2;
        uniDexBest123 = [uniDexBest12];
    else
        uniDexBest3 = unique(ShallowDailyTempBest3(:,1));
        uniDexBest123 = [uniDexBest12; uniDexBest3];
    end
    if length(uniDexBest123) == length(uniDex)
        ShallowDailyTempBest = [ShallowDailyTempBest1; ShallowDailyTempBest2; ShallowDailyTempBest3];
        save -mat DailyShallow_83052_Best.mat ShallowDailyTempBest DepthTagMax
    end
    % elseif length(uniDexBest12) == length(uniDex)
    %     ShallowDailyTempBest = ShallowDailyTempBest1;
    %     save -mat DailyShallow_83052_Best.mat ShallowDailyTempBest DepthTagMax
end

%% Fourth pass: 25

a = exist('uniDexBest123', 'var'); % since most won't get to the fourth pass
if a == 1
    if length(uniDexBest123) < length(uniDex)
        toCompleteDex4 = []; % figure out which of uniDex are still to be completed
        for i = 1:length(uniDex)
            ilog = uniDexBest123 == uniDex(i);
            if sum(ilog) <= 0 % if this is not in uniDexBest4
                toCompleteDex4 = [toCompleteDex4; uniDex(i)];
            end
        end
        
        load -mat DailyShallow_83052_25.mat
        co = 25;
        
        uniDexCount4 = [];
        ShallowDailyTempBest4 = [];
        for i = 1:length(toCompleteDex4)
            ilog = toCompleteDex4(i) == ShallowDailyTemp(:,1);
            uniDexCount4 = [uniDexCount4 sum(ilog)];
            if sum(ilog) >= nMin
                ShallowDailyTempBest4 = [ShallowDailyTempBest4; ShallowDailyTemp(ilog,:) repmat(co,sum(ilog),1)];
            end
        end
        a = isempty(ShallowDailyTempBest4);
        if a == 1 % troubleshoot: if the 15m increment doesn't add anything
            ShallowDailyTempBest4 = ShallowDailyTempBest3;
            uniDexBest1234 = uniDexBest123;
        else
            uniDexBest4 = unique(ShallowDailyTempBest4(:,1));
            uniDexBest1234 = [uniDexBest123; uniDexBest4];
        end
        if length(uniDexBest1234) == length(uniDex)
            ShallowDailyTempBest = [ShallowDailyTempBest1; ShallowDailyTempBest2; ShallowDailyTempBest3; ShallowDailyTempBest4];
            save -mat DailyShallow_83052_Best.mat ShallowDailyTempBest DepthTagMax
        end
        % elseif length(uniDexBest123) == length(uniDex)
        %     ShallowDailyTempBest = ShallowDailyTempBest1;
        %     save -mat DailyShallow_83052_Best.mat ShallowDailyTempBest DepthTagMax
    end
    
    %% Fourth pass: 45
    a = exist('uniDexBest1234', 'var'); % since most won't get to the fourth pass
    if a == 1
        if length(uniDexBest1234) < length(uniDex)
            toCompleteDex5 = []; % figure out which of uniDex are still to be completed
            for i = 1:length(uniDex)
                ilog = uniDexBest1234 == uniDex(i);
                if sum(ilog) <= 0 % if this is not in uniDexBest5
                    toCompleteDex5 = [toCompleteDex5; uniDex(i)];
                end
            end
            
            load -mat DailyShallow_83052_50.mat % will change this to 45 for tag 83050+1 and +2, % 60 for 83040+6
            co = 50;
            
            uniDexCount5 = [];
            ShallowDailyTempBest5 = [];
            for i = 1:length(toCompleteDex5)
                ilog = toCompleteDex5(i) == ShallowDailyTemp(:,1);
                uniDexCount5 = [uniDexCount5 sum(ilog)];
                if sum(ilog) >= nMin
                    ShallowDailyTempBest5 = [ShallowDailyTempBest5; ShallowDailyTemp(ilog,:) repmat(co,sum(ilog),1)];
                end
            end
            a = isempty(ShallowDailyTempBest5);
            if a == 1 % troubleshoot: if the 15m increment doesn't add anything
                ShallowDailyTempBest5 = ShallowDailyTempBest4;
                uniDexBest12345 = uniDexBest1234;
            else
                uniDexBest5 = unique(ShallowDailyTempBest5(:,1));
                uniDexBest12345 = [uniDexBest1234; uniDexBest5];
            end
            if length(uniDexBest12345) == length(uniDex)
                ShallowDailyTempBest = [ShallowDailyTempBest1; ShallowDailyTempBest2; ShallowDailyTempBest3; ShallowDailyTempBest4; ShallowDailyTempBest5];
                save -mat DailyShallow_83052_Best.mat ShallowDailyTempBest DepthTagMax
            end
        end
    end
end

% and one final troubleshoot
a = exist('ShallowDailyTempBest', 'var');
if a == 0
    if tagNum == 83050+2 || tagNum == 83040+6 % write it like this so search and replace doesn' muck it up
        if length(uniDex)-length(uniDexBest12345) == 1
            toCompleteDex6 = []; % figure out which of uniDex are still to be completed
            for i = 1:length(uniDex)
                ilog = uniDexBest12345 == uniDex(i);
                if sum(ilog) <= 0 % if this is not in uniDexBest5
                    toCompleteDex6 = uniDex(i);
                end
            end
            ilog2 = toCompleteDex6 == ShallowDailyTemp(:,1);
            ShallowDailyTempBest6 = [ShallowDailyTemp(ilog2,:) repmat(co,sum(ilog2),1)];
            ShallowDailyTempBest = [ShallowDailyTempBest1; ShallowDailyTempBest2; ShallowDailyTempBest3;...
                ShallowDailyTempBest4; ShallowDailyTempBest5; ShallowDailyTempBest6];
            save -mat DailyShallow_83052_Best.mat ShallowDailyTempBest DepthTagMax
        end
    else
        error('do another shallow depth cutoff');
    end
end
%%


%%

disp('Completed determineShallowCutoff.m')
% ===== EOF [determineShallowCutoff.m] ======
