% rd_MSReboundVBehavBin.m

%% settings
experiment = 'e0e3e5'; % 'e0','e3','e5','e0e3e5'
analysisName = 'bonneh_post2_dir'; % 'bonneh_pre2','bonneh_post2','bonneh_pre2_dir','bonneh_post2_dir'

fileName = sprintf('data/%s/LatencyVBehavAnalysis_%s.mat', experiment, analysisName);

saveFileName = sprintf('data/%s/LatencyVBehavBinAnalysis_%s.mat', experiment, analysisName);
saveAnalysis = 0;

normalizeOption = 'none'; % 'none','DC','morey'

if strfind(analysisName, 'dir')
    analyzeDir = 1;
else
    analyzeDir = 0;
end

%% load data
load(fileName)

%% setup
if strfind(analysisName,'pre')
    win = [500 1500]; % pre2
elseif strfind(analysisName,'post')
    win = [1500 2700]; % post2
else
    error('analysisName not recognized')
end

binSize = 200; %200, 100
binStep = 200; %200, 10
% win = win + 100; %win-binSize/2; % shift if needed, + for pre, - for post
% binStarts = win(1):binSize:win(2)-1;
% binStarts = win(1):binSize/2:win(2)-1-binSize/2; % overlapping bins
binStarts = win(1):binStep:win(2)-binSize; % overlapping bins
nBins = numel(binStarts);

nDirBins = 16; % direction bins

if strcmp(experiment,'e5')
    targetNames = {'T1','T2','T3'};
    targetTimes = [1500 1750 2000];
    cueNames = {'n','t1','t2','t3'};
else
    targetNames = {'T1','T2'};
    targetTimes = [1500 1750];
    cueNames = {'n','t1','t2'};
end
nTargets = numel(targetNames);
nCues = numel(cueNames);

nSubjects = numel(subjects);

%% get data
% make everything the same size
sz = size(msData);
accData = accData(1:sz(1),:,:);
targetCond = targetCond(1:sz(1),:,:);
cueCond = repmat(1:3,sz(1),1,sz(3));
cueCond(isnan(accData)) = NaN;

% forget about cue condition (long format)
msDataL = reshape(msData,sz(1)*sz(2),sz(3));
accDataL = reshape(accData,sz(1)*sz(2),sz(3));
targetCondL = reshape(targetCond,sz(1)*sz(2),sz(3));
cueCondL = reshape(cueCond,sz(1)*sz(2),sz(3));
if analyzeDir
    dirDataL = reshape(dirData,sz(1)*sz(2),sz(3));
end

%% calculate accuracy and ms direction as a function of bin and target
nTrials = [];
for iBin = 1:nBins
    bin = [binStarts(iBin) binStarts(iBin)+binSize];
    wBin = msDataL>=bin(1) & msDataL<bin(2);
    
    for iT = 1:nTargets
        wT = targetCondL==iT;
        
        w = wBin & wT;
        nTrials(iBin,iT,:) = sum(w,1);
        
        for iS = 1:nSubjects
            accBin(iBin,iT,iS) = nanmean(accDataL(w(:,iS),iS),1);
            if analyzeDir
                dirBin{iBin,iT,iS} = dirDataL(w(:,iS),iS);
                dirT{iT,iS} = dirDataL(wT(:,iS),iS);
            end
        end
    end
end

%% calculate accuracy and ms direction as a function of bin and cue
nTrials = [];
for iBin = 1:nBins
    bin = [binStarts(iBin) binStarts(iBin)+binSize];
    wBin = msDataL>=bin(1) & msDataL<bin(2);
    
    for iC = 1:nCues
        wC = cueCondL==iC;
        
        w = wBin & wC;
        nTrialsCue(iBin,iC,:) = sum(w,1);
        
        for iS = 1:nSubjects
            accBinCue(iBin,iC,iS) = nanmean(accDataL(w(:,iS),iS),1);
            if analyzeDir
                dirBinCue{iBin,iC,iS} = dirDataL(w(:,iS),iS);
                dirC{iC,iS} = dirDataL(wC(:,iS),iS);
            end
        end
    end
end

%% summary stats
accBinMean = nanmean(accBin,3); % nan ! for post2

% normalized error bars
switch normalizeOption
    case 'none'
        accBinSte = nanstd(accBin,0,3)/sqrt(nSubjects);
    case 'DC'
        accBinNorm = normalizeDC(accBin);
        accBinSte = nanstd(accBinNorm,0,3)/sqrt(nSubjects);
    case 'morey'
        m = 2; % M/(M-1), use 2 bc comparing each bin to mean (essentially seprate t-tests)
        accBinNorm = normalizeDC(accBin); 
        accBinSte = nanstd(accBinNorm,0,3)/sqrt(nSubjects)*sqrt(m); 
    otherwise
        error('normalizeOption not recognized')
end

%% calculate accuracy as a function of target only
nTrialsT = [];
for iT = 1:nTargets
    w = targetCondL==iT;
    nTrialsT(iT,:) = sum(w,1);
    
    for iS = 1:nSubjects
        accT(iT,iS) = mean(accDataL(w(:,iS),iS),1);
    end
end

accTMean = mean(accT,2);
accTSte = std(accT,0,2)./sqrt(nSubjects);

%% bin accuracy difference from mean accuracy
accTBin(1,:,:) = accT;
accTBin = repmat(accTBin,nBins,1,1);

accBinDiff = accBin - accTBin;
accBinDiffMean = nanmean(accBinDiff,3);
accBinDiffSte = nanstd(accBinDiff,0,3)./sqrt(nSubjects);

%% stats
[h p] = ttest(accBin(:,1,:),accBin(:,2,:),'dim',3)

hM = []; pM = [];
for iT = 1:nTargets
    [hM(:,iT) pM(:,iT)] = ttest(squeeze(accBin(:,iT,:)),...
        repmat(accT(iT,:),nBins,1),'dim',2);
end

%% ms direction distribution as a function of bin and target
if analyzeDir
% dirBin{iBin,iT,iS}
% direction setup
dirBinSize = 2*pi/nDirBins;
dirBinEdges = (-dirBinSize:dirBinSize:2*pi-dirBinSize) + dirBinSize/2;

% count directions
counts = zeros(length(dirBinEdges)-1, nBins, nTargets, nSubjects);
rates = zeros(size(counts));
props = zeros(size(counts));

countsT = zeros(length(dirBinEdges)-1, nTargets, nSubjects);
ratesT = zeros(size(countsT));

for iS = 1:nSubjects
    for iT = 1:nTargets
        nTrialsInTargetCond = nnz(targetCondL(:,iS)==iT);
        dT = dirT{iT,iS};
        hT = polarhistogram(dT,dirBinEdges);
        countsT(:,iT,iS) = hT.Values;
        ratesT(:,iT,iS) = hT.Values/sum(nTrials(:,iT,iS));
        for iBin = 1:nBins
            d = dirBin{iBin,iT,iS};
            h = polarhistogram(d,dirBinEdges);
            counts(:,iBin,iT,iS) = h.Values;
            rates(:,iBin,iT,iS) = h.Values/nTrials(iBin,iT,iS); % /trial % (binSize/1000*nTrials(iBin,iT,iS))
            props(:,iBin,iT,iS) = h.Values/nTrialsInTargetCond;
        end
    end
end

% average counts and rates per direction
groupDirCountMean = mean(counts,4);
groupDirRateMean = mean(rates,4);
groupDirPropMean = mean(props,4);

groupDirCountMeanT = mean(countsT,3);
groupDirRateMeanT = mean(ratesT,3);

groupDirRateMeanAll = mean(squeeze(sum(countsT,2))./repmat(squeeze(sum(sum(nTrials,1),2))',nDirBins,1),2);
end

%% ms direction distribution as a function of bin and cue
if analyzeDir
% dirBinCue{iBin,iC,iS}
% direction setup
dirBinSize = 2*pi/nDirBins;
dirBinEdges = (-dirBinSize:dirBinSize:2*pi-dirBinSize) + dirBinSize/2;

% count directions
countsCue = zeros(length(dirBinEdges)-1, nBins, nCues, nSubjects);
ratesCue = zeros(size(countsCue));
propsCue = zeros(size(countsCue));

countsC = zeros(length(dirBinEdges)-1, nCues, nSubjects);
ratesC = zeros(size(countsC));

for iS = 1:nSubjects
    for iC = 1:nCues
        nTrialsInCueCond = nnz(cueCondL(:,iS)==iC);
        dC = dirC{iC,iS};
        hC = polarhistogram(dC,dirBinEdges);
        countsC(:,iC,iS) = hC.Values;
        ratesC(:,iC,iS) = hC.Values/sum(nTrialsCue(:,iC,iS));
        for iBin = 1:nBins
            d = dirBinCue{iBin,iC,iS};
            h = polarhistogram(d,dirBinEdges);
            countsCue(:,iBin,iC,iS) = h.Values;
            ratesCue(:,iBin,iC,iS) = h.Values/nTrialsCue(iBin,iC,iS); % /trial % (binSize/1000*nTrials(iBin,iT,iS))
            propsCue(:,iBin,iC,iS) = h.Values/nTrialsInCueCond;
        end
    end
end

% average counts and rates per direction
groupDirCountCueMean = mean(countsCue,4);
groupDirRateCueMean = nanmean(ratesCue,4);
groupDirPropCueMean = mean(propsCue,4);

groupDirCountMeanC = mean(countsC,3);
groupDirRateMeanC = mean(ratesC,3);

groupDirRateCueMeanAll = mean(squeeze(sum(countsC,2))./repmat(squeeze(sum(sum(nTrialsCue,1),2))',nDirBins,1),2);
end


%% plot setup
figure
colors = get(gca,'ColorOrder');
close(gcf)

% % colors(3,:) = [.5 .5 .5];
% c = rgb2hsv(colors);
% c(:,2) = .2; % saturation
% c(:,3) = .95; % value
% % c(3,3) = .9;
% ebcolors = hsv2rgb(c);

% colors = [colors(4,:); .98 .56 .02];
ebcolors = [colors repmat(.4,size(colors,1),1)]; % add alpha column
grayBoxColor = [.9 .9 .9];

targetDur = 30;

%% plot
figure
histogram(msDataL(:))
xlabel('latency (ms)')
ylabel('number of trials')

figure
hold on
for iS = 1:nSubjects
    plot(repmat(binStarts',1,nTargets)+binSize/2, nTrials(:,:,iS))
    set(gca,'ColorOrderIndex',1) % reset color order
end
xlabel('bin center (ms)')
ylabel('number of trials')
legend(targetNames)

figure
hold on
for iS = 1:nSubjects
    plot(repmat(binStarts',1,nTargets)+binSize/2, accBin(:,:,iS))
    set(gca,'ColorOrderIndex',1) % reset color order
%     fprintf('%d\n',iS)
%     pause(1)
end
xlabel('bin center (ms)')
ylabel('proportion correct')
legend(targetNames)

xlims = [binStarts(1) binStarts(end)+binSize];
ylims = [.45 .8];
if nTargets==3
    j = [-10 0 10];
else
    j = [-10 10];
end
figure
hold on
for iT = 1:nTargets
    plot(xlims,[accTMean(iT) accTMean(iT)],'Color',ebcolors(iT,:),'LineWidth',4)
end
for iT = 1:nTargets
    errbar(binStarts+binSize/2+j(iT),accBinMean(:,iT),accBinSte(:,iT),...
        'Color',colors(iT,:))
    p1(iT) = plot(binStarts+binSize/2+j(iT),accBinMean(:,iT),...
        'Color',colors(iT,:),'Marker','.','MarkerSize',30);
end
xlabel('Latency bin (ms)')
ylabel('Proportion correct')
title(und2space(analysisName))
xlim(xlims)
ylim(ylims)
legend(p1, targetNames)
legend boxoff
set(gca,'TickDir','out','LineWidth',1,'XTick',binStarts+binSize/2)

% diff from mean
xlims = [binStarts(1) binStarts(end)+binSize];
ylims = [-.1501 .1501];
if nTargets==3
    j = [-10 0 10];
else
    j = [-10 10];
end
figure
hold on
plotboxes([targetTimes' (targetTimes+repmat(targetDur,1,nTargets))'],ylims, grayBoxColor)
plot(xlims,[0 0],'k')
for iT = 1:nTargets
    errbar(binStarts+binSize/2+j(iT),accBinDiffMean(:,iT),accBinDiffSte(:,iT),...
        'Color',colors(iT,:))
    p1(iT) = plot(binStarts+binSize/2+j(iT),accBinDiffMean(:,iT),...
        'Color',colors(iT,:),'Marker','.','MarkerSize',30);
end
xlabel('Latency bin (ms)')
ylabel('Proportion correct, difference from mean')
title(und2space(analysisName))
xlim(xlims)
ylim(ylims)
legend(p1, targetNames)
legend boxoff
set(gca,'TickDir','out','LineWidth',1,'XTick',binStarts+binSize/2)

% diff from mean - line with shaded error bars 
xlims = [binStarts(1) binStarts(end)+binSize];
% ylims = [-.1501 .1501];
ylims = [-.25 .25];
figure
hold on
% plotboxes([targetTimes' (targetTimes+repmat(targetDur,1,nTargets))'],ylims, grayBoxColor)
plotboxes([1560 1710], ylims, grayBoxColor)
for iEvent = 1:numel(targetTimes)
    plot([targetTimes(iEvent) targetTimes(iEvent)],ylims,'--k')
end
plot(xlims,[0 0],'k')
for iT = 1:nTargets
    shadedErrorBar(binStarts+binSize/2,accBinDiffMean(:,iT),accBinDiffSte(:,iT),...
        {'Color',colors(iT,:)},0)
    p1(iT) = plot(binStarts+binSize/2,accBinDiffMean(:,iT),...
        'Color',colors(iT,:));
end
xlabel('Latency bin (ms)')
ylabel('Proportion correct, difference from mean')
title(und2space(analysisName))
xlim(xlims)
ylim(ylims)
legend(p1, targetNames)
legend boxoff
set(gca,'TickDir','out','LineWidth',1)
% pbaspect([750/500 1 1])

%% plot directions
if analyzeDir
% every direction, regardless of subject (so some subjects will count more)
figure
polarhistogram(dirDataL(:),'BinEdges',dirBinEdges)

% mean across subjects, all conditions
groupDirCountMeanAll = sum(sum(groupDirCountMean,3),2);
figure
polarhistogram('BinEdges',dirBinEdges,'BinCounts',groupDirRateMeanAll)
title('proportion of trials')
set(gca,'FontSize',16)
set(gca,'ThetaTick',0:22.5:360)
set(gca,'RLim',[0 .2])

% mean per condition, regardless of time bin
vals = groupDirRateMeanT;
rlims = [0 .3]; %[0 70];

figure
for iT = 1:nTargets
    h = polarhistogram('BinEdges',dirBinEdges,'BinCounts',vals(:,iT))';
    h.DisplayStyle = 'stairs';
    h.LineWidth = 2;
    h.Parent.RLim = rlims;
    hold on
end
title('proportion of trials')
legend(targetNames)

% mean per condition
measure = 'props'; % 'counts', 'rates', 'props'
switch measure
    case 'counts'
        vals = groupDirCountMean;
        rlims = [0 20]; 
        titleText = 'number of MS';
    case 'rates'
        vals = groupDirRateMean;
        rlims = [0 .4]; 
        titleText = 'proportion of trials';
    case 'props'
        vals = groupDirPropMean;
        rlims = [0 .03]; %[0 .07];
        titleText = 'p(MS)';
    otherwise
        error('measure not recognized')
end

figure('Position',[150 650 2000 550])
for iBin = 1:nBins
    subplot(1,nBins+1,iBin)
    for iT = 1:nTargets
        h = polarhistogram('BinEdges',dirBinEdges,'BinCounts',vals(:,iBin,iT));
        h.DisplayStyle = 'stairs';
        h.LineWidth = 2;
        h.Parent.RLim = rlims;
        hold on
    end
    title(sprintf('%d-%d ms',binStarts(iBin),binStarts(iBin)+binSize))
end
% plot legend
subplot(1,nBins+1,iBin+1)
for iT = 1:nTargets
    h = polarhistogram('BinEdges',dirBinEdges,'BinCounts',ones(1,length(dirBinEdges)-1));
    h.DisplayStyle = 'stairs';
    h.LineWidth = 2;
    h.Parent.RLim = rlims;
    hold on
end
title(titleText)
legend(targetNames)
end

%% %% plot directions - by cue
if analyzeDir
% mean per condition
measure = 'props'; % 'counts', 'rates', 'props'
switch measure
    case 'counts'
        vals = groupDirCountCueMean;
        rlims = [0 20]; 
        titleText = 'number of MS';
    case 'rates'
        vals = groupDirRateCueMean;
        rlims = [0 .4]; 
        titleText = 'proportion of trials';
    case 'props'
        vals = groupDirPropCueMean;
        rlims = [0 .04]; 
        titleText = 'p(MS)';
    otherwise
        error('measure not recognized')
end

figure('Position',[150 650 2000 550])
for iBin = 1:nBins
    subplot(1,nBins+1,iBin)
    for iC = 1:nCues
        h = polarhistogram('BinEdges',dirBinEdges,'BinCounts',vals(:,iBin,iC));
        h.DisplayStyle = 'stairs';
        h.LineWidth = 2;
        h.Parent.RLim = rlims;
        hold on
    end
    title(sprintf('%d-%d ms',binStarts(iBin),binStarts(iBin)+binSize))
end
% plot legend
subplot(1,nBins+1,iBin+1)
for iC = 1:nCues
    h = polarhistogram('BinEdges',dirBinEdges,'BinCounts',ones(1,length(dirBinEdges)-1));
    h.DisplayStyle = 'stairs';
    h.LineWidth = 2;
    h.Parent.RLim = rlims;
    hold on
end
title(titleText)
legend(cueNames)
end

%% save
if saveAnalysis
    save(saveFileName, 'subjects', 'exp', 'condNames', 'accs', 'accNames', ...
        'targets', 'targetNames', 'normalizeOption', ...
        'win', 'binStarts', 'binSize', 'nBins', 'targetTimes', ...
        'nTrials', 'accBin', 'accBinMean', 'accBinSte', ...
        'nTrialsT', 'accT', 'accTMean', 'accTSte', 'accTBin', ...
        'accBinDiff', 'accBinDiffMean', 'accBinDiffSte')
end
