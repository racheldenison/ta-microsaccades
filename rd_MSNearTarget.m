% rd_MSNearTarget.m

experiment = 'e5';

saveFileName = sprintf('data/%s/LatencyNearTargetAnalysis.mat', experiment);
saveAnalysis = 0;

targetWin = [-100 0]; 

switch experiment
    case {'e0','e3','e5'}
        load(sprintf('data/%s/onsets.mat',experiment))
        load(sprintf('data/%s/behav_data.mat',experiment))
        
    case 'e0e3e5'
        e0Onsets = load('data/e0/onsets.mat');
        e3Onsets = load('data/e3/onsets.mat');
        e5Onsets = load('data/e5/onsets.mat');
        condsAll = cat(2, e0Onsets.condsAll, e3Onsets.condsAll, e5Onsets.condsAll);
        onsetsAll = cat(2, e0Onsets.onsetsAll, e3Onsets.onsetsAll, e5Onsets.onsetsAll);
        
        e0Behav = load('data/e0/behav_data.mat');
        e3Behav = load('data/e3/behav_data.mat');
        e5Behav = load('data/e5/behav_data.mat');
        
        sz = size(e5Behav.accData);
        sz1 = size(e0Behav.accData);
        accData1 = nan(sz(1),sz(2),sz1(3));
        accData1(1:sz1(1),1:sz1(2),1:end) = e0Behav.accData;
        rcc1 = nan(sz(1),sz(2),sz1(3));
        rcc1(1:sz1(1),1:sz1(2),1:end) = e0Behav.respCueCond;
        
        sz2 = size(e3Behav.accData);
        accData2 = nan(sz(1),sz(2),sz2(3));
        accData2(1:sz2(1),1:sz2(2),1:end) = e3Behav.accData;
        rcc2 = nan(sz(1),sz(2),sz2(3));
        rcc2(1:sz1(1),1:sz1(2),1:end) = e3Behav.respCueCond;
        
        accData3 = e5Behav.accData;
        rcc3 = e5Behav.respCueCond;
        
        accData = cat(3, accData1, accData2, accData3);
        respCueCond = cat(3, rcc1, rcc2, rcc3);
        
        condNames = e0Onsets.condNames;
        
    otherwise
        error('experiment not recognized')
end

nSubjects = numel(onsetsAll);
t = 1:size(onsetsAll{1},2);

switch experiment
    case {'e0','e3','e0e3e5'}
        targetNames = {'t1','t2'};
        targetTimes = [1500 1750];
        cueNames = {'n','t1','t2'};
    case 'e5'
        targetNames = {'t1','t2','t3'};
        targetTimes = [1500 1750 2000];
        cueNames = {'n','t1','t2','t3'};
    otherwise
        error('experiment not recognized')
end

nCues = numel(cueNames);
nTargets = numel(targetNames);

%% reorganize accData and targetCond
for iS = 1:nSubjects
    conds = condsAll{iS};
    nConds = numel(unique(conds));
    for iCond = 1:nConds
        nTrialsPerCueCond(iCond) = nnz(conds==iCond);
    end
    accDataAll{iS} = [];
    targetCondAll{iS} = [];
    for iCond = 1:nConds
        accDataAll{iS} = [accDataAll{iS}; accData(1:nTrialsPerCueCond(iCond),iCond,iS)];
        targetCondAll{iS} = [targetCondAll{iS}; respCueCond(1:nTrialsPerCueCond(iCond),iCond,iS)];
    end
end

%% find trials with onsets within each target window
for iS = 1:nSubjects
    onsets = onsetsAll{iS};

    for iT = 1:nTargets
        onsetsNearTarget = onsets(:,targetTimes(iT)+targetWin(1):targetTimes(iT)+targetWin(2));
        hasOnsetNearTarget{iS}(:,iT) = any(onsetsNearTarget,2);
        nOnsetNearTarget{iS}(:,iT) = sum(onsetsNearTarget,2);
    end
end

%% acc depending on whether a MS was near the target
hs = [0 1];
accNT = [];
for iS = 1:nSubjects
    target = targetCondAll{iS};
    acc = accDataAll{iS};
    hasONT = hasOnsetNearTarget{iS};
    
    for iT = 1:nTargets
        wT = target==iT;
        
        for iH = 1:numel(hs)
            h = logical(hs(iH));
            wH = hasONT(:,iT)==h;
            w = wT & wH;
            nTrialsNT(iH,iT,iS) = nnz(w);
            
            accNT(iH,iT,iS) = nanmean(acc(w));
        end
    end
end

accNTMean = nanmean(accNT,3);
accNTSte = nanstd(accNT,0,3)./sqrt(nSubjects);

accNTDiff = squeeze(accNT(2,:,:) - accNT(1,:,:));
accNTDiffMean = nanmean(accNTDiff,2);
accNTDiffSte = nanstd(accNTDiff,0,2)./sqrt(nSubjects);

disp(mean(nTrialsNT,3))

figure
bar(squeeze(nTrialsNT(2,:,:))')
xlabel('subject')
ylabel('number of trials')
legend(targetNames)

figure
hold on
colors = get(gca,'ColorOrder');
for iT = 1:nTargets
    p1(iT,:) = plot(squeeze(accNT(:,iT,:)),'.-','Color',colors(iT,:),'MarkerSize',20);
end
legend(p1(:,1),targetNames)
xlabel('MS near target')
ylabel('proportion correct')

figure
errorbar(repmat([0; 1],1,nTargets), accNTMean, accNTSte)
set(gca,'XTick',hs)
legend(targetNames)
xlabel('MS near target')
ylabel('proportion correct')

figure
hold on
plot([.8 nTargets+.2],[0 0],'--k')
errorbar(1:nTargets, accNTDiffMean, accNTDiffSte,'.','MarkerSize',20)
set(gca,'XTick',1:nTargets,'XTickLabel',targetNames)
xlabel('target')
ylabel('\Delta proportion correct')
xlim([.8 nTargets+.2])
pbaspect([.5 1 1])

%% number of MS near target as a function of whether the target was cued
msNT0 = [];
for iS = 1:nSubjects
    cue = condsAll{iS};
    nONT = nOnsetNearTarget{iS};
    
    for iC = 1:nCues
        w = cue==iC;
        for iT = 1:nTargets
            msNTCount(iC,iT,iS) = sum(nONT(w,iT));
            msNT0(iC,iT,iS) = mean(nONT(w,iT));
        end
    end
end

% grab a subset
msNT = [];
for iT = 1:nTargets
    msNT(1,iT,:) = msNT0(strcmp(cueNames,'n'),iT,:); % cue neutral
    msNT(2,iT,:) = msNT0(strcmp(cueNames,targetNames{iT}),iT,:); % cue target
end

msNTMean = nanmean(msNT,3);
msNTSte = nanstd(msNT,0,3)./sqrt(nSubjects);

figure
errorbar(msNTMean, msNTSte)
set(gca,'XTick',1:nCues,'XTickLabel',{'neutral','cue target'})
legend(targetNames)
ylabel('proportion trials with MS')

figure
errorbar(msNTMean', msNTSte','.','MarkerSize',20)
set(gca,'XTick',1:nTargets,'XTickLabel',targetNames)
xlabel('target')
ylabel('proportion trials with MS')
legend({'neutral','cue target'},'box','off')
xlim([.8 nTargets+.2])
pbaspect([.5 1 1])
box off
jitterx(gca);

%% save
if saveAnalysis
    save(saveFileName, 'experiment', 'targetWin', 'condNames', ...
        'accData', 'respCueCond', 't', 'targetNames', 'targetTimes', 'cueNames', ...
        'accDataAll', 'targetCondAll','hasOnsetNearTarget','nOnsetNearTarget',...
        'hs','accNT','nTrialsNT','accNTMean','accNTSte',... 
        'accNTDiff','accNTDiffMean','accNTDiffSte',...
        'msNTCount', 'msNT0', 'msNT', 'msNTMean', 'msNTSte')
end
