% rd_MSReboundVBehav.m

%% settings
experiment = 'e5'; % 'e0','e3','e5','e0e3e5'
analysisName = 'bonneh_pre2_dir'; %'bonneh','bonneh_pre','bonneh_pre2','bonneh_post','bonneh_post2'

saveFileName = sprintf('data/%s/LatencyVBehavAnalysis_%s.mat', experiment, analysisName);
saveAnalysis = 0;

accNames = {'correct','incorrect'};
accs = [1 0]; 
nAccs = numel(accs);

if strcmp(experiment,'e5')
    condNames = {'n','t1','t2','t3'};
else
    condNames = {'n','t1','t2'};
end
nConds = numel(condNames);

targets = 1:nConds-1;
nTargets = numel(targets);
targetNames = condNames(2:end);

%% load data
switch experiment
    case {'e0','e3','e5'}
        behavFile = sprintf('data/%s/behav_data.mat', experiment);
        b = load(behavFile);
        accData = b.accData;
        targetCond = b.respCueCond;
        
        msFile = sprintf('data/%s/%s.mat', experiment, analysisName);
        ms = load(msFile);
        msData = ms.lastOnset;
        dirData = ms.lastDir;
 
        load(sprintf('data/%s/subjects.mat', experiment))
        exp = repmat(cellstr(experiment),1,numel(subjects));
        
    case 'e0e3e5'
        e0B = load('data/e0/behav_data.mat');
        e3B = load('data/e3/behav_data.mat');
        e5B = load('data/e5/behav_data.mat');
        
        e0MS = load(sprintf('data/e0/%s.mat', analysisName));
        e3MS = load(sprintf('data/e3/%s.mat', analysisName));
        e5MS = load(sprintf('data/e5/%s.mat', analysisName));
        
        e0S = load('data/e0/subjects.mat');
        e3S = load('data/e3/subjects.mat');
        e5S = load('data/e5/subjects.mat');
        subjects = [e0S.subjects e3S.subjects e5S.subjects];
        exp = [repmat(cellstr('e0'),1,numel(e0S.subjects)) ...
            repmat(cellstr('e3'),1,numel(e3S.subjects)) ...
            repmat(cellstr('e5'),1,numel(e5S.subjects))];
        
        sz = max([size(e0B.accData,1), size(e3B.accData,1), size(e5B.accData,1)]);
        accData = nan(sz, size(e0B.accData,2), numel(subjects));
        accData(1:size(e0B.accData,1),:,strcmp(exp,'e0')) = e0B.accData;
        accData(1:size(e3B.accData,1),:,strcmp(exp,'e3')) = e3B.accData;
        accData(1:size(e5B.accData,1),:,strcmp(exp,'e5')) = e5B.accData(:,1:3,:);
        
        sz = max([size(e0B.respCueCond,1), size(e3B.respCueCond,1), size(e5B.respCueCond,1)]);
        targetCond = nan(sz, size(e0B.respCueCond,2), numel(subjects));
        targetCond(1:size(e0B.respCueCond,1),:,strcmp(exp,'e0')) = e0B.respCueCond;
        targetCond(1:size(e3B.respCueCond,1),:,strcmp(exp,'e3')) = e3B.respCueCond;
        targetCond(1:size(e5B.respCueCond,1),:,strcmp(exp,'e5')) = e5B.respCueCond(:,1:3,:);
        
        sz = max([size(e0MS.lastOnset,1), size(e3MS.lastOnset,1), size(e5MS.lastOnset,1)]);
        msData = nan(sz, size(e0MS.lastOnset,2), numel(subjects));
        msData(1:size(e0MS.lastOnset,1),:,strcmp(exp,'e0')) = e0MS.lastOnset;
        msData(1:size(e3MS.lastOnset,1),:,strcmp(exp,'e3')) = e3MS.lastOnset;
        msData(1:size(e5MS.lastOnset,1),:,strcmp(exp,'e5')) = e5MS.lastOnset(:,1:3,:);
        
        sz = max([size(e0MS.lastDir,1), size(e3MS.lastDir,1), size(e5MS.lastDir,1)]);
        dirData = nan(sz, size(e0MS.lastDir,2), numel(subjects));
        dirData(1:size(e0MS.lastDir,1),:,strcmp(exp,'e0')) = e0MS.lastDir;
        dirData(1:size(e3MS.lastDir,1),:,strcmp(exp,'e3')) = e3MS.lastDir;
        dirData(1:size(e5MS.lastDir,1),:,strcmp(exp,'e5')) = e5MS.lastDir(:,1:3,:);
        
    otherwise
        error('experiment not recognized')
end

nSubjects = numel(subjects);

%% summary stats
groupTrialData = [];
groupData = [];
groupTrialDataT = [];
groupDataT = [];
nTrialsT = [];
nTrialsTWithData = [];

for iS = 1:nSubjects
    for iCond = 1:nConds
        accVals = accData(:,iCond,iS);
        targetVals = targetCond(:,iCond,iS);
        msVals = msData(:,iCond,iS);
        
        for iAcc = 1:nAccs
            acc = accs(iAcc);
            
            wAcc = accVals==acc;
            groupTrialData{iCond,iAcc,iS} = msVals(wAcc);
            groupData(iCond,iAcc,iS) = nanmean(msVals(wAcc));
            
            for iT = 1:nTargets
                target = targets(iT);
                
                wT = targetVals==target;
                w = wAcc & wT;
                
                groupTrialDataT{iCond,iAcc,iT,iS} = msVals(w);
                groupDataT(iCond,iAcc,iT,iS) = nanmean(msVals(w));
                nTrialsT(iCond,iAcc,iT,iS) = nnz(w);
                nTrialsTWithData(iCond,iAcc,iT,iS) = nnz(~isnan(msVals(w)));
            end
        end
    end
end

groupMean = mean(groupData,3);
groupSte = std(groupData,0,3)./sqrt(nSubjects);

groupMeanT = nanmean(groupDataT,4); % nan!
groupSteT = nanstd(groupDataT,0,4)./sqrt(nSubjects);

%% by acc, regardless of condition
for iS = 1:nSubjects
    trialDataByAcc{1,iS} = [];
    trialDataByAcc{2,iS} = [];

    for iCond = 1:nConds
        for iAcc = 1:nAccs
            trialDataByAcc{iAcc,iS} = [trialDataByAcc{iAcc,iS}; groupTrialData{iCond,iAcc,iS}];
            groupDataByAcc(iAcc,iS) = nanmean(trialDataByAcc{iAcc,iS});
        end
    end
end

groupMeanByAcc = mean(groupDataByAcc,2);
groupSteByAcc = std(groupDataByAcc,0,2)./sqrt(nSubjects);

%% by acc and response cue, average across cue condition
groupDataTCueAve = squeeze(nanmean(groupDataT,1)); % more nan!
groupMeanTCueAve = mean(groupDataTCueAve,3);
groupSteTCueAve = std(groupDataTCueAve,0,3)./sqrt(nSubjects);

%% plot
xlims = [.8 nConds+.2];

%% number of trials
figure
histogram(nTrialsT(:),'BinMethod','integer')
xlabel('number of trials in condition')
ylabel('number of conditions, across subjects')

figure
histogram(nTrialsTWithData(:),'BinMethod','integer')
xlabel('number of trials in condition contributing a microsaccade')
ylabel('number of conditions, across subjects')

%% by acc and cue condition
figure
for iSubject = 1:nSubjects
    subplot(nSubjects,1,iSubject)
    hold on
    plot(groupData(:,:,iSubject))
    xlim(xlims)
    set(gca,'XTick',1:nConds)
    set(gca,'XTickLabel',[])
    if iSubject==1
        legend(accNames)
    elseif iSubject==nSubjects
        set(gca,'XTickLabel',condNames)
        xlabel('cue')
        ylabel('latency (ms)')
    end
    yyaxis right
    set(gca,'YTick',[])
    ylabel(subjects{iSubject})
end

figure
errorbar(groupMean, groupSte)
legend(accNames)
set(gca,'XTick',1:nConds)
set(gca,'XTickLabel',condNames)
xlim(xlims)
xlabel('cue')
ylabel('latency (ms)')
title(und2space(analysisName))

%% by acc, cue condition, and response cue condition
lineStyles = {'-','--',':'};
colors = get(gca,'ColorOrder');
idx = 1; names = [];
figure
hold on
for iAcc = 1:nAccs
    for iT = 1:nTargets
        errorbar(groupMeanT(:,iAcc,iT), groupSteT(:,iAcc,iT),'Color',colors(iAcc,:),'LineStyle',lineStyles{iT})
        names{idx} = sprintf('target %s, %s', condNames{iT+1}, accNames{iAcc});
        idx = idx + 1;
    end
end
set(gca,'XTick',1:nConds)
set(gca,'XTickLabel',condNames)
xlim(xlims)
xlabel('cue')
ylabel('latency (ms)')
legend(names)
title(und2space(analysisName))

% good figure
figure
hold on
for iT = 1:nTargets
    subplot(1,nTargets,iT)
    errorbar(groupMeanT(:,:,iT), groupSteT(:,:,iT))
    set(gca,'XTick',1:nConds)
    set(gca,'XTickLabel',condNames)
    xlim(xlims)
    xlabel('cue')
    ylabel('latency (ms)')
    title(targetNames{iT})
end
legend(accNames)
rd_supertitle2(und2space(analysisName))

lineStyles = {'-','--',':','-.'};
colors = get(gca,'ColorOrder');
idx = 1; names = [];
figure
hold on
for iCond = 1:nConds
    for iAcc = 1:nAccs
        errorbar(squeeze(groupMeanT(iCond,iAcc,:)), squeeze(groupSteT(iCond,iAcc,:)),'Color',colors(iAcc,:),'LineStyle',lineStyles{iCond})
        names{idx} = sprintf('cue %s, %s', condNames{iCond}, accNames{iAcc});
        idx = idx + 1;
    end
end
set(gca,'XTick',1:nTargets)
set(gca,'XTickLabel',condNames(2:end))
xlim(xlims)
xlabel('target')
ylabel('latency (ms)')
legend(names)
title(und2space(analysisName))

%% by acc and response cue, average across condition
figure
errorbar(groupMeanTCueAve', groupSteTCueAve')
set(gca,'XTick',1:nTargets)
set(gca,'XTickLabel',condNames(2:end))
xlabel('target')
ylabel('latency (ms)')
xlim([.8 2.2])
legend(accNames)
title(und2space(analysisName))

ylims = [min(groupMeanTCueAve(:))*.995 max(groupMeanTCueAve(:))*1.005];
figure
bar(groupMeanTCueAve')
ylim(ylims)
set(gca,'XTickLabel',condNames(2:end))
xlabel('target')
ylabel('latency (ms)')
legend(accNames)
title(und2space(analysisName))

%% by acc, regardless of condition
figure
plot(groupDataByAcc)
set(gca,'XTick',1:nAccs)
set(gca,'XTickLabel',accNames)
ylabel('latency (ms)')
title(und2space(analysisName))

figure
errorbar(groupMeanByAcc, groupSteByAcc)
set(gca,'XTick',1:nAccs)
set(gca,'XTickLabel',accNames)
ylabel('latency (ms)')
title(und2space(analysisName))

%% save
groupDataTCueAve = squeeze(nanmean(groupDataT,1)); % more nan!
groupMeanTCueAve = mean(groupDataTCueAve,3);
groupSteTCueAve = std(groupDataTCueAve,0,3)./sqrt(nSubjects);
if saveAnalysis
    save(saveFileName, 'subjects', 'exp', 'condNames', 'accs', 'accNames', 'targets', ...
        'accData', 'targetCond', 'msData', 'dirData', ...
        'groupTrialData', 'groupData', 'groupMean', 'groupSte', ...
        'groupTrialDataT', 'groupDataT', 'groupMeanT', 'groupSteT', ...
        'nTrialsT', 'nTrialsTWithData', ...
        'trialDataByAcc', 'groupDataByAcc', 'groupMeanByAcc', 'groupSteByAcc', ...
        'groupDataTCueAve', 'groupMeanTCueAve', 'groupSteTCueAve')
end

