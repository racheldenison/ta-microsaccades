% rd_MSRebound.m

%% settings
experiment = 'e0e3e5'; %'e0e3e5_nothresh'; %'e0','e3','e5','e0e3e5','e0_nothresh',etc,'e0e3e5_nothresh'
analysisName = 'bonneh_pre2'; %'bonneh','bonneh_pre','bonneh_pre2','bonneh_post','bonneh_post2'
removeSubjects = [];
removeRepeatSubjects = false; 

saveFileName = sprintf('data/%s/LatencyAnalysis_%s.mat', experiment, analysisName);
saveAnalysis = 0;

%% load data
switch experiment
    case {'e0','e3','e5','e0_nothresh','e3_nothresh','e5_nothresh'}
        dataFile = sprintf('data/%s/%s.mat', experiment, analysisName);
        load(dataFile)
        data = boneh_mat;
        if exist('lastOnset','var')
            trialsData = lastOnset;
        else
            trialsData = [];
        end
        load(sprintf('data/%s/subjects.mat', experiment))
        exp = repmat(cellstr(experiment),1,numel(subjects));
        
    case 'e0e3e5'
        e0 = load(sprintf('data/e0/%s.mat', analysisName));
        e3 = load(sprintf('data/e3/%s.mat', analysisName));
        e5 = load(sprintf('data/e5/%s.mat', analysisName));
        data = cat(1, e0.boneh_mat, e3.boneh_mat, e5.boneh_mat(:,1:3));
        
        e0S = load('data/e0/subjects.mat');
        e3S = load('data/e3/subjects.mat');
        e5S = load('data/e5/subjects.mat');
        subjects = [e0S.subjects e3S.subjects e5S.subjects];
        exp = [repmat(cellstr('e0'),1,numel(e0S.subjects)) ...
            repmat(cellstr('e3'),1,numel(e3S.subjects)) ...
            repmat(cellstr('e5'),1,numel(e5S.subjects))];
        
        sz = max([size(e0.lastOnset,1), size(e3.lastOnset,1), size(e5.lastOnset,1)]);
        trialsData = nan(sz, size(e0.lastOnset,2), numel(subjects));
        trialsData(1:size(e0.lastOnset,1),:,strcmp(exp,'e0')) = e0.lastOnset;
        trialsData(1:size(e3.lastOnset,1),:,strcmp(exp,'e3')) = e3.lastOnset;
        trialsData(1:size(e5.lastOnset,1),:,strcmp(exp,'e5')) = e5.lastOnset(:,1:3,:);
        
    case 'e0e3e5_nothresh'
        e0 = load(sprintf('data/e0_nothresh/%s.mat', analysisName));
        e3 = load(sprintf('data/e3_nothresh/%s.mat', analysisName));
        e5 = load(sprintf('data/e5_nothresh/%s.mat', analysisName));
        data = cat(1, e0.boneh_mat, e3.boneh_mat, e5.boneh_mat(:,1:3));
        
        e0S = load('data/e0_nothresh/subjects.mat');
        e3S = load('data/e3_nothresh/subjects.mat');
        e5S = load('data/e5_nothresh/subjects.mat');
        subjects = [e0S.subjects e3S.subjects e5S.subjects];
        exp = [repmat(cellstr('e0'),1,numel(e0S.subjects)) ...
            repmat(cellstr('e3'),1,numel(e3S.subjects)) ...
            repmat(cellstr('e5'),1,numel(e5S.subjects))];
        
        n = max([size(e0.lastOnset,1), size(e3.lastOnset,1), size(e5.lastOnset,1)]);
        trialsData = nan(n, size(e0.lastOnset,2), numel(subjects));
        trialsData(1:size(e0.lastOnset,1),:,strcmp(exp,'e0')) = e0.lastOnset;
        trialsData(1:size(e3.lastOnset,1),:,strcmp(exp,'e3')) = e3.lastOnset;
        trialsData(1:size(e5.lastOnset,1),:,strcmp(exp,'e5')) = e5.lastOnset(:,1:3,:);
    otherwise
        error('experiment not recognized')
end

%% optionally remove subjects
rsidx = [];
if ~isempty(removeSubjects)
    for iS = 1:numel(removeSubjects)
        idx = find(strcmp(subjects, removeSubjects{iS}));
        if ~isempty(idx)
            rsidx(iS) = idx;
        end
    end
    rsidx(rsidx==0) = [];
end

if removeRepeatSubjects
    [us, usidx] = unique(subjects);
    rsidx = setdiff(1:numel(subjects),usidx);
end

%% setup
switch experiment
    case {'e0','e3','e0e3e5','e0_nothresh','e3_nothresh','e0e3e5_nothresh'}
        condNames = {'n','t1','t2'};
    case {'e5','e5_nothresh'}
        condNames = {'n','t1','t2','t3'};
    otherwise
        error('experiment not found')
end

nSubjects = size(data,1);
nConds = numel(condNames);

data(rsidx,:) = NaN;

%% summary stats
groupMean = nanmean(data,1);
groupSte = nanstd(data,0,1)./sqrt(nSubjects-numel(rsidx));

dataNorm = normalizeDC(data')';
groupMeanNorm = nanmean(dataNorm,1);
groupSteNorm = nanstd(dataNorm,0,1)./sqrt(nSubjects-numel(rsidx));

%% trial median and summary stats
dataM = squeeze(nanmedian(trialsData,1))';

groupMMean = nanmean(dataM,1);
groupMSte = nanstd(dataM,0,1)./sqrt(nSubjects-numel(rsidx));

dataMNorm = normalizeDC(dataM')';
groupMMeanNorm = nanmean(dataMNorm,1);
groupMSteNorm = nanstd(dataMNorm,0,1)./sqrt(nSubjects-numel(rsidx));

for iS = 1:nSubjects
    vals = trialsData(:,:,iS);
    trialsDataNoCond(:,iS) = vals(:);
end
dataMNoCond = nanmedian(trialsDataNoCond,1);
[min(dataMNoCond-1500) max(dataMNoCond-1500)]

nT1DiffM = groupMMean(strcmp(condNames,'n'))-groupMMean(strcmp(condNames,'t1'));

%% normalized trialsData
for iSubject = 1:nSubjects
    vals = trialsData(:,:,iSubject);
    valsX = nanmean(vals(:));
    valsSigma = nanstd(vals(:));
    trialsDataXShift(:,:,iSubject) = vals - valsX;
    trialsDataZScore(:,:,iSubject) = (vals - valsX)./valsSigma;
end

xgridz = linspace(-5,5,100);
trialsDataZScoreDensity = [];
for iSubject = 1:nSubjects
    for iCond = 1:nConds
        vals = trialsDataZScore(:,iCond,iSubject);
        [trialsDataZScoreDensity(:,iCond,iSubject)] = ksdensity(vals, xgridz);
    end
end

% sample equal number of trials per subject and condition
rng(1,'twister'); % fix random number generator
n = 100; % number of sampled trials
for iSubject = 1:nSubjects
    for iCond = 1:nConds
        vals = trialsDataZScore(:,iCond,iSubject);
        valsn = vals(~isnan(vals));
        idx = randi(numel(valsn),n,1); % bootstrap sample
        trialsDataZScoreSampled(:,iCond,iSubject) = valsn(idx);
    end
end

%% summary stats on normalized data
dataZ = squeeze(nanmedian(trialsDataZScore,1))';
groupZMean = nanmean(dataZ,1);
groupZSte = nanstd(dataZ,0,1)./sqrt(nSubjects-numel(rsidx));

%% F test
[fvals, pF, rn, tbl] = rd_rmANOVA(data, condNames, {'Attention'}, nConds);
disp(tbl)

%% t-test, KS test
switch experiment
    case {'e0','e3','e0e3e5','e0_nothresh','e3_nothresh','e0e3e5_nothresh'}
        for iCond = 1:nConds
            if iCond==nConds
                pairCond = 1;
            else
                pairCond = iCond+1;
            end
            str = sprintf('%s%s', condNames{iCond}, condNames{pairCond});
            vals = squeeze(data(:,iCond)-data(:,pairCond))';
            [hT.(str), pT.(str), ~, stat] = ttest(vals);
            tvals.(str) = stat.tstat;
            
            vals = squeeze(dataM(:,iCond)-dataM(:,pairCond))';
            [hTM.(str), pTM.(str), ~, stat] = ttest(vals);
            tvalsM.(str) = stat.tstat;
            
            vals = squeeze(dataZ(:,iCond)-dataZ(:,pairCond))';
            [hTZ.(str), pTZ.(str), ~, stat] = ttest(vals);
            tvalsZ.(str) = stat.tstat;
            
            vals1 = squeeze(trialsDataZScoreSampled(:,iCond,:));
            vals2 = squeeze(trialsDataZScoreSampled(:,pairCond,:));
            [hKS.(str), pKS.(str), statKS.(str)] = kstest2(vals1(:), vals2(:));
        end
    case {'e5','e5_nothresh'}
        for iCond = 2:nConds
            if iCond==nConds
                pairCond = 2;
            else
                pairCond = iCond+1;
            end
            str = sprintf('%s%s', condNames{iCond}, condNames{pairCond});
            vals = squeeze(data(:,iCond)-data(:,pairCond))';
            [hT.(str), pT.(str), ~, stat] = ttest(vals);
            tvals.(str) = stat.tstat;
            
            vals = squeeze(dataM(:,iCond)-dataM(:,pairCond))';
            [hTM.(str), pTM.(str), ~, stat] = ttest(vals);
            tvalsM.(str) = stat.tstat;
            
            vals = squeeze(dataZ(:,iCond)-dataZ(:,pairCond))';
            [hTZ.(str), pTZ.(str), ~, stat] = ttest(vals);
            tvalsZ.(str) = stat.tstat;
            
            vals1 = squeeze(trialsDataZScoreSampled(:,iCond,:));
            vals2 = squeeze(trialsDataZScoreSampled(:,pairCond,:));
            [hKS.(str), pKS.(str), statKS.(str)] = kstest2(vals1(:), vals2(:));
        end
        for iCond = 2:nConds
            pairCond = 1;
            str = sprintf('%s%s', condNames{iCond}, condNames{pairCond});
            vals = squeeze(data(:,iCond)-data(:,pairCond))';
            [hT.(str), pT.(str), ~, stat] = ttest(vals);
            tvals.(str) = stat.tstat;
            
            vals = squeeze(dataM(:,iCond)-dataM(:,pairCond))';
            [hTM.(str), pTM.(str), ~, stat] = ttest(vals);
            tvalsM.(str) = stat.tstat;
            
            vals = squeeze(dataZ(:,iCond)-dataZ(:,pairCond))';
            [hTZ.(str), pTZ.(str), ~, stat] = ttest(vals);
            tvalsZ.(str) = stat.tstat;
            
            vals1 = squeeze(trialsDataZScoreSampled(:,iCond,:));
            vals2 = squeeze(trialsDataZScoreSampled(:,pairCond,:));
            [hKS.(str), pKS.(str), statKS.(str)] = kstest2(vals1(:), vals2(:));
        end
end

%% count ms in certain time windows
refTime = 1750;
for iS = 1:nSubjects
    vals = trialsData(:,:,iS);
    trialsDataBySubject(:,iS) = vals(:);
    msProp(iS) = nnz(trialsDataBySubject(:,iS)<refTime)/nnz(~isnan(trialsDataBySubject(:,iS)));
end
msPropAll = nnz(trialsDataBySubject(:)<refTime)/nnz(~isnan(trialsDataBySubject(:)));

%% plot
figure
bar(data)
legend(condNames)
xlabel('observer')
ylabel('rebound latency (ms)')

figure
errorbar(groupMean, groupSte, '.', 'MarkerSize',20)
xlim([0.5 numel(condNames)+0.5])
set(gca,'XTick',1:numel(condNames))
set(gca,'XTickLabel',condNames)
ylabel('rebound latency (ms)')

figure
errorbar(groupMeanNorm, groupSteNorm, '.', 'MarkerSize',20)
xlim([0.5 numel(condNames)+0.5])
set(gca,'XTick',1:numel(condNames))
set(gca,'XTickLabel',condNames)
ylabel('rebound latency (ms)')

figure
bar(dataM)
legend(condNames)
xlabel('observer')
ylabel('rebound latency trial median (ms)')

figure
errorbar(groupMMean, groupMSte, '.', 'MarkerSize',20)
xlim([0.5 numel(condNames)+0.5])
set(gca,'XTick',1:numel(condNames))
set(gca,'XTickLabel',condNames)
ylabel('rebound latency trial median (ms)')

figure
errorbar(groupMMeanNorm, groupMSteNorm, '.', 'MarkerSize',20)
xlim([0.5 numel(condNames)+0.5])
set(gca,'XTick',1:numel(condNames))
set(gca,'XTickLabel',condNames)
ylabel('rebound latency trial median (ms)')

figure
bar(dataZ)
legend(condNames)
xlabel('observer')
ylabel('rebound latency (z-score)')

figure
errorbar(groupZMean, groupZSte, '.', 'MarkerSize',20)
xlim([0.5 numel(condNames)+0.5])
set(gca,'XTick',1:numel(condNames))
set(gca,'XTickLabel',condNames)
ylabel('rebound latency (z-score)')

figure
plot([data(:,strcmp(condNames,'t1')), data(:,strcmp(condNames,'t2'))]')

if nConds==3
    pairs = {'t1','t2'; 't1','n'; 'n','t2'};
else
    pairs = {'t1','t2'; 't1','n'; 'n','t2'; 't1','t3'; 'n','t3'};
end
nPairs = size(pairs,1);
exps = unique(exp);
lims = [min(data(:))/1.1 max(data(:))*1.1];
figure('Position',[100 300 950 300])
colors = get(gca,'ColorOrder');
for iPair = 1:nPairs
    pair = pairs(iPair,:);
    subplot(1,nPairs,iPair)
    hold on
    for iExp = 1:numel(exps)
        inExp = strcmp(exp, exps{iExp});
        plot(data(inExp,strcmp(condNames,pair{1})), ...
            data(inExp,strcmp(condNames,pair{2})),'o','Color',colors(iExp,:))
    end
    plot(lims,lims,'k')
    axis equal
    xlim(lims)
    ylim(lims)
    xlabel(pair(1))
    ylabel(pair(2))
end
legend(exps,'Location','best')

if numel(exps)==3
    expNums = [ones(1,numel(e0S.subjects)) ...
        2*ones(1,numel(e3S.subjects)) ...
        3*ones(1,numel(e5S.subjects))];
    if strfind(analysisName, 'pre')
        ylims = [700 1000];
    else
        ylims = [1500 3000];
    end
    figure
    for iExp=1:3
        subplot(1,3,iExp)
        bar(mean(data(expNums==iExp,:)));
        ylim(ylims)
        set(gca,'XTickLabel',condNames)
        title(exps{iExp})
    end
end

if ~isempty(trialsData)
    normalizeOption = 'zscoresamp';
    switch normalizeOption
        case 'xshift'
            tData = trialsDataXShift;
            xlimsPre = [-800 800];
            binWidth = 50;
        case 'zscore'
            tData = trialsDataZScore;
            xlimsPre = [-5 5];
            xlimsPost = [-5 5];
            binWidth = 0.25;
        case 'zscoresamp'
            tData = trialsDataZScoreSampled;
            xlimsPre = [-5 5];
            xlimsPost = [-5 5];
            binWidth = 0.25;
        case 'none'
            tData = trialsData;
            xlimsPre = [500 1500];
%             xlimsPost = [1500 4500];
            xlimsPost = [1500 3000];
            binWidth = 50;
        otherwise
            error('normalizeOption not recognized')
    end
    if strfind(analysisName, 'pre')
        xlims = xlimsPre;
    else
        xlims = xlimsPost;
    end
    
    figure
    for iSubject = 1:nSubjects
        [h(iSubject),p(iSubject),ks2stat(iSubject)] = ...
            kstest2(trialsData(:,1,iSubject), trialsData(:,2,iSubject));
        subplot(nSubjects,1,iSubject)
        hold on
        for iCond = 1:nConds
            histogram(tData(:,iCond,iSubject),'BinWidth',binWidth,'Normalization','pdf')
        end
        if iSubject==1
            legend(condNames)
        elseif iSubject==nSubjects
            xlabel('latency (ms index)')
            ylabel('probability density')
        end
        xlim(xlims)
        yyaxis right
        set(gca,'YTick',[])
        ylabel(subjects{iSubject})
    end
    
    figure
    subplot(2,1,1)
    hold on
    for iCond = 1:nConds
        vals = tData(:,iCond,:);
        histogram(vals(:),'BinWidth',binWidth,'Normalization','pdf')
    end
    xlabel('latency (ms index)')
    ylabel('probability density')
    legend(condNames)
    subplot(2,1,2)
    hold on
    xgrid = linspace(xlims(1),xlims(2),100);
    tDataDensity = [];
    for iCond = 1:nConds
        vals = tData(:,iCond,:);
        [tDataDensity(:,iCond)] = ksdensity(vals(:), xgrid);
        plot(xgrid,tDataDensity(:,iCond),'color',colors(iCond,:))
    end
    xlabel('latency (ms index)')
    ylabel('probability density')
    legend(condNames)
    
    figure
    hold on
    for iCond = 1:nConds
        shadedErrorBar(xgridz, mean(trialsDataZScoreDensity(:,iCond,:),3), ...
            std(trialsDataZScoreDensity(:,iCond,:),0,3)./sqrt(nSubjects), ...
            {'Color',colors(iCond,:)},1);
    end
    xlim([-3 3])
    xlabel('latency (z-score)')
    ylabel('probability density')
    
    msInTrial = ~isnan(trialsData);
    msInTrialCount = squeeze(sum(msInTrial,1));
    figure
    bar(msInTrialCount')
    xlabel('subject')
    ylabel('number of trials with microsaccade in time period')
    legend(condNames)
end

figure
bar(msProp)
xlabel('subject')
ylabel(sprintf('proportion MS < %d ms', refTime))

%% save
if saveAnalysis
    save(saveFileName, 'subjects', 'condNames', 'data', 'groupMean', 'groupSte', ...
        'dataNorm', 'groupMeanNorm', 'groupSteNorm', ...
        'dataM', 'groupMMean', 'groupMSte', 'dataMNorm', 'groupMMeanNorm', 'groupMSteNorm', ...
        'dataZ', 'groupZMean', 'groupZSte', ...
        'hT', 'pT', 'tvals', 'hTM', 'pTM', 'tvalsM', 'hTZ', 'pTZ', 'tvalsZ', ...
        'trialsData', 'trialsDataXShift', 'trialsDataZScore', 'trialsDataZScoreDensity', ...
        'xgridz', 'msInTrial', 'msInTrialCount')
end

