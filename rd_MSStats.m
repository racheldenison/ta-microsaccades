% rd_MSStats.m

%% settings
experiment = 'e0e3e5'; %'e0','e3','e5','e0e3','e0e3e5','e0_nothresh',etc,'e0e3e5_nothresh'
removeSubjects = [];
% removeSubjects = {'rd','jp','ec'};
removeRepeatSubjects = false;

cueColors = {
    [.5 .5 .5]
    [123 146 202]./255
    [242 123 96]./255
    [.49 .18 .55]};

yellowShadeColor = [1 .93 .8];

saveFileName = sprintf('data/%s/RateAnalysis.mat', experiment);
saveAnalysis = false;

%% load data
switch experiment
    case 'e0'
        load data/e0/analysis_struct4.mat
        load(sprintf('data/%s/subjects.mat', experiment))
    case 'e3'
        load data/e3/analysis_struct-11_subjects.mat
        load(sprintf('data/%s/subjects.mat', experiment)) 
    case 'e5'
        load('data/e5/analysis_struct.mat');
        load(sprintf('data/%s/subjects.mat', experiment)) 
    case 'e0e3'
        e0 = load('data/e0/analysis_struct4.mat');
        e3 = load('data/e3/analysis_struct-11_subjects.mat'); 
        analysis_struct.results_per_subject = ...
            [e0.analysis_struct.results_per_subject, e3.analysis_struct.results_per_subject];
        e0S = load(sprintf('data/%s/subjects.mat', 'e0')); 
        e3S = load(sprintf('data/%s/subjects.mat', 'e3')); 
        subjects = [e0S.subjects e3S.subjects];
    case 'e0e3e5'
        e0 = load('data/e0/analysis_struct4.mat');
        e3 = load('data/e3/analysis_struct-11_subjects.mat');
        e5 = load('data/e5/analysis_struct.mat');
        analysis_struct.results_per_subject = ...
            [e0.analysis_struct.results_per_subject, e3.analysis_struct.results_per_subject, ...
            e5.analysis_struct.results_per_subject];
        e0S = load(sprintf('data/%s/subjects.mat', 'e0'));
        e3S = load(sprintf('data/%s/subjects.mat', 'e3'));
        e5S = load(sprintf('data/%s/subjects.mat', 'e5'));
        subjects = [e0S.subjects e3S.subjects e5S.subjects];
    case 'e0e3e5_nothresh'
        e0 = load('data/e0_nothresh/analysis_struct.mat');
        e3 = load('data/e3_nothresh/analysis_struct.mat');
        e5 = load('data/e5_nothresh/analysis_struct.mat');
        analysis_struct.results_per_subject = ...
            [e0.analysis_struct.results_per_subject, e3.analysis_struct.results_per_subject, ...
            e5.analysis_struct.results_per_subject];
        e0S = load(sprintf('data/%s/subjects.mat', 'e0_nothresh'));
        e3S = load(sprintf('data/%s/subjects.mat', 'e3_nothresh'));
        e5S = load(sprintf('data/%s/subjects.mat', 'e5_nothresh'));
        subjects = [e0S.subjects e3S.subjects e5S.subjects];
    otherwise
        load(sprintf('data/%s/analysis_struct.mat', experiment))
        load(sprintf('data/%s/subjects.mat', experiment))
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
measure = 'microsaccadic_rate';

dataStruct = analysis_struct.results_per_subject;
nSubjects = numel(dataStruct);

conds = fields(dataStruct{1}.(measure));
nConds = numel(conds);

% cut off end of time series
switch experiment
    case {'e5','e0e3e5','e5_nothresh','e0e3e5_nothresh'}
        nSamplesAllowed = 2700;
        for iSubject = 1:nSubjects
            for iCond = 1:nConds
                dataStruct{iSubject}.(measure).(conds{iCond})(:,nSamplesAllowed+1:end) = [];
            end
        end
end

nSamples = length(dataStruct{1}.(measure).(conds{1}));

% The time window for the first 2 experiments was [-1500, 2250], and for 
% the third it was [-1500, 3000]. I created though the microsaccades rate 
% graphs with an averaging window with a width of 50ms, so it cut off 25ms 
% from the start and 25ms from the end of the data. --Omer
switch experiment
    case {'e0','e3','e0e3','e0e3e5','e0_nothresh','e3_nothresh','e0e3e5_nothresh'}
%         t = -1499:1200; % check
%         t = -1474:2225; % check
        t = -1474:1225; % check
        condNames = {'n','t1','t2'};
        eventTimes = [0 250];
    case {'e5','e5_nothresh'}
%         t = -1474:2975;
        t = -1474:1225;
        condNames = {'n','t1','t2','t3'};
        eventTimes = [0 250 500];
    otherwise
        error('experiment not found')
end

toi = -499:0; % -750:750
toiIdx = find(t==toi(1)):find(t==toi(end));
nT = numel(toi);

% check data dimention
if nSamples ~= numel(t)
    error('check number of samples')
end
if ~isempty(subjects) && nSubjects ~= numel(subjects)
    error('check number of subjects')
end

%% organize data
data = zeros(nSamples,nConds,nSubjects);

for iSubject = 1:nSubjects
    for iCond = 1:nConds
        data(:,iCond,iSubject) = dataStruct{iSubject}.(measure).(conds{iCond});
    end
end
data(:,:,rsidx) = NaN;

% summary
dataMean = nanmean(data,3);
dataSte = nanstd(data,0,3)./sqrt(nSubjects-numel(rsidx)); 

%% load stats
% clusters.betaCueNeutral = dlmread('data/combined/TimelinesPermutationTests_msrate_beta_cueneutral_noheader.txt');
% clusters.betaCueT2 = dlmread('data/combined/TimelinesPermutationTests_msrate_beta_cuet2_noheader.txt');
if strfind(experiment, 'nothresh')
    clusters.betaCueNeutral = dlmread('data/combined/TimelinesPermutationTests_E0E3E5_msrate_-500to0ms_N30_nothresh_beta_cueneutral_noheader.txt');
    clusters.betaCueT2 = [];
else
    clusters.betaCueNeutral = dlmread('data/combined/TimelinesPermutationTests_E0E3E5_msrate_-500to0ms_N30_beta_cueneutral_noheader.txt');
    clusters.betaCueT2 = [];
end

%% basic time series plots
figure('Position',[150 800 1100 450])
hold on
plot(t, dataMean)
ylim([0 3])
for i = 1:nConds-1
    vline(eventTimes(i),'Color','k','LineStyle',':','LineWidth',1.5)
end
legend(condNames)
xlabel('time (ms)')
ylabel('microsaccade rate (Hz)')

figure('Position',[500 150 750 1100])
for iSubject = 1:nSubjects
    subplot(nSubjects,1,iSubject)
    hold on
    plot(t, data(:,:,iSubject))
    xlim([t(1) t(end)])
    ylim([0 5])
    for i = 1:nConds-1
        vline(eventTimes(i),'Color','k','LineStyle',':','LineWidth',1.5)
    end
    yyaxis right
    set(gca,'YTick',[])
    ylabel(subjects{iSubject})
    if iSubject==1
        legend(condNames)
    end
    if iSubject~=nSubjects
        set(gca,'XTick',[])
    end
end
xlabel('time (ms)')
yyaxis left
ylabel('microsaccade rate (Hz)')

%% time series with shaded significant clusters
clusterNames = fields(clusters);
clusterColors = {[.8 .8 .8], [251 197 185]./255};
plotClusters = 1;

ylims = [0 1.8];
figure('Position',[150 800 1100 450])
hold on
% rectangle('Position',[-500, ylims(1), 500, diff(ylims)], 'FaceColor', [.95 .95 .95], 'EdgeColor', 'none')
if plotClusters
    for iComp = 1:numel(clusterNames)
        clusterName = clusterNames{iComp};
        vals = clusters.(clusterName);
        %     vals = vals(vals(:,1)==max(vals(:,1)),:);
        for iC = 1:size(vals,1)
            clusterTimes = vals(iC,2:3);
            rectangle('Position',[clusterTimes(1), ylims(1), diff(clusterTimes), diff(ylims)],...
                'FaceColor', clusterColors{iComp}, 'EdgeColor', 'none')
        end
    end
end
for i = 1:nConds
    plot(t, dataMean(:,i), 'Color', cueColors{i},'LineWidth',6)
end
ylim(ylims)
for i = 1:nConds-1
    vline(eventTimes(i),'Color','k','LineStyle',':','LineWidth',2.5)
end
legend(condNames)
xlim([-500 0])
xlabel('Time (ms)')
ylabel('Microsaccade rate (Hz)')

legend('neutral','cue T1','cue T2')
legend boxoff
set(gca,'LineWidth',1)
set(gca,'TickDir','out')
set(gca,'FontSize',20)

% print_pdf('msrate_rebound', '~/Desktop')

%% get data from clusters
clusterName = 'betaCueNeutral';
clusterVals = clusters.(clusterName);
nClusters = size(clusterVals,1);

clusterData = [];
for iC = 1:nClusters
    clusterTimes = clusterVals(iC,2:3);
    tidx = find(t==clusterTimes(1)):find(t==clusterTimes(2));
    clusterData(:,iC,:) = mean(data(tidx,:,:),1);
end

clusterMean = mean(clusterData,3);
clusterSte = std(clusterData,0,3)./sqrt(nSubjects);

% for e5
% for iC = 1:nClusters
%     [h p] = ttest(clusterData(2,iC,:), clusterData(4,iC,:))
% end

figure
errorbar(clusterMean, clusterSte)
xlim([.5 nConds+.5])
set(gca,'XTick',1:nConds)
set(gca,'XTickLabel',condNames)
ylabel('Microsaccade rate (Hz)')

%% check for unequal variance
% vartestn(x) returns a summary table of statistics and a box plot for a 
% Bartlett test of the null hypothesis that the columns of data vector x 
% come from normal distributions with the same variance.
for iSample = 1:nSamples
    vals = squeeze(data(iSample,:,:))';
    pBartlett(iSample) = vartestn(vals,'Display','off');
    pOBrien(iSample) = vartestn(vals,'Display','off','TestType','OBrien');
end

figure
hold on
plot(t, pBartlett)
plot(t, pOBrien)
plot([t(1) t(end)],[.05 .05],'k')
xlabel('time (ms)')
ylabel('p-value')
legend('Bartlett','OBrien')

%% check for normally distributed data
for iSample = 1:nSamples
    vals = squeeze(data(iSample,:,:))';
    hLillie(iSample) = lillietest(vals(:));
    hJB(iSample) = jbtest(vals(:));
end

figure
hold on
plot(t, hLillie)
plot(t, hJB + 0.1)
xlabel('time (ms)')
legend('Lillie','JB')

%% empirical F-test
alpha = 0.05;
tic
for iSample = 1:nSamples
    vals = squeeze(data(iSample,:,:))';
    [fvals(iSample), pF(iSample)] = rd_rmANOVA(vals, conds, {'Attention'}, nConds);
end
toc

[fvalsClusterSums, fvalsMaxAbsClusterSum, fvalsClusterWins] = rd_clusterSum(fvals(toiIdx), pF(toiIdx)<alpha);

% find max cluster window
fvalsMaxClusterWin = ...
    toi(fvalsClusterWins(abs(fvalsClusterSums)==fvalsMaxAbsClusterSum,:));


figure
plot(t, fvals)
xlabel('time (ms)')
ylabel('F stat')

figure
hold on
plot(t, pF)
plot(t, pBartlett)
plot([t(1) t(end)],[alpha alpha],'k')
ylim([0 2*alpha])
xlabel('time (ms)')
ylabel('p-value')
legend('F test','Bartlett test')

figure
hold on
plot(t, pF<alpha)
plot(t, (pOBrien<alpha) + 0.1)
ylim([-1 2])
xlabel('time (ms)')
ylabel(sprintf('p-value < %.2f', alpha))
legend('F test','OBrien test')
set(gca,'YTickLabel',[])

%% empirical Friedman's test (nonparametric 1-way repeated-measures ANOVA)
tic
for iSample = 1:nSamples
    vals = squeeze(data(iSample,:,:))';
    [pFriedman(iSample),tbl] = friedman(vals,1,'off');
    chisqFriedman(iSample) = tbl{2,5};
end
toc

[friedmanClusterSums, friedmanMaxAbsClusterSum] = rd_clusterSum(chisqFriedman(toiIdx), pFriedman(toiIdx)<alpha);

figure
hold on
plot(t, fvals)
plot(t, chisqFriedman)
xlabel('time (ms)')
legend('F stat','Friedman chi-sq')

figure
hold on
plot(t, pF)
plot(t, pFriedman)
plot([t(1) t(end)],[.05 .05],'k')
xlabel('time (ms)')
legend('F stat','Friedman chi-sq')

figure
hold on
plot(t, pF<alpha)
plot(t, (pFriedman<alpha) + 0.1)
ylim([-1 2])
xlabel('time (ms)')
ylabel(sprintf('p-value < %.2f', alpha))
legend('F test','Friedman test')
set(gca,'YTickLabel',[])

%%
% [p,tbl,stats] = kruskalwallis(x,[],'off');
% [p,tbl,stats] = friedman(x,reps,displayopt);

%% t-test
alpha = 0.05;
for iCond = 1:nConds %2:nConds
    if iCond==nConds
        pairCond = 1; %2
    else
        pairCond = iCond+1; %1
    end
    str = sprintf('%s%s', condNames{iCond}, condNames{pairCond});
    vals = squeeze(data(:,iCond,:)-data(:,pairCond,:))';
    [hT.(str), pT.(str), ~, stat] = ttest(vals);
    tvals.(str) = stat.tstat;
end

m = fieldnames(hT);
for iM = 1:numel(m)
    [tvalsClusterSums.(m{iM}), tvalsMaxAbsClusterSum.(m{iM}), tvalsClusterWins.(m{iM})] = ...
        rd_clusterSum(tvals.(m{iM})(toiIdx), pT.(m{iM})(toiIdx)<alpha);
end

% find max cluster window
for iM = 1:numel(m)
    tvalsMaxClusterWin.(m{iM}) = ...
        toi(tvalsClusterWins.(m{iM})(abs(tvalsClusterSums.(m{iM}))==tvalsMaxAbsClusterSum.(m{iM}),:));
end

figure
subplot(2,1,1)
hold on
for iM = 1:numel(m)
    plot(t, hT.(m{iM})+0.1*(iM-1))
end
ylabel('t stat h')
legend(m)
subplot(2,1,2)
hold on
for iM = 1:numel(m)
    plot(t, pT.(m{iM}))
end
xlabel('time (ms)')
ylabel('t stat p')

%% permutation F-stat null distribution
nShuffles = 1000;

fvalsSh = []; pvalsSh = [];
fvalsClusterSumsSh = []; fvalsMaxAbsClusterSumSh = [];
tic
for iSh = 1:nShuffles
    fprintf('\nShuffle %d\n', iSh)
    for iSubject = 1:nSubjects
       dataSh(:,:,iSubject) = data(:,randperm(nConds),iSubject); 
    end
    parfor iT = 1:nT
        if mod(iT,100)==0, fprintf('.'), end
        if mod(iT,1000)==0, fprintf('\n'), end
        vals = squeeze(dataSh(toiIdx(iT),:,:))';
        [fvalsSh(iT,iSh), pvalsSh(iT,iSh), rowNames] = ...
            rd_rmANOVA(vals, conds, {'Attention'}, nConds);
    end
    
    [fvalsClusterSumsSh{iSh}, fvalsMaxAbsClusterSumSh(iSh)] = ...
        rd_clusterSum(fvalsSh(:,iSh), pvalsSh(:,iSh)<alpha);
    
    if mod(iSh,10)==0, toc, end
end

% confidence intervals
fvalsCI = prctile(fvalsSh,[2.5 97.5],2);

% cluster sum p-value
clusterPVal = nnz(fvalsMaxAbsClusterSumSh>fvalsMaxAbsClusterSum)/numel(fvalsMaxAbsClusterSumSh);

% plot
figure
hold on
plot(t, fvals)
plot(toi, fvalsCI, 'k')
ylabel('F')
title('permutation F-stat null distribution')

figure
hold on
histogram(fvalsMaxAbsClusterSumSh)
vline(fvalsMaxAbsClusterSum,'k')
ax = gca;
ht = text(fvalsMaxAbsClusterSum*1.1,ax.YLim(2)*0.9,sprintf('p = %.2f', clusterPVal));
ht.FontSize = 16;
xlabel('F stat max cluster sum')

%% permutation Friedman null distribution
nShuffles = 1000;

friedmanSh = []; pFriedmanSh = [];
friedmanClusterSumsSh = []; friedmanMaxAbsClusterSumSh = [];
tic
for iSh = 1:nShuffles
    fprintf('\nShuffle %d\n', iSh)
    for iSubject = 1:nSubjects
       dataSh(:,:,iSubject) = data(:,randperm(nConds),iSubject); 
    end
    parfor iT = 1:nT
        if mod(iT,100)==0, fprintf('.'), end
        if mod(iT,1000)==0, fprintf('\n'), end
        vals = squeeze(dataSh(toiIdx(iT),:,:))';
        [pFriedmanSh(iT,iSh),tbl] = friedman(vals,1,'off');
        friedmanSh(iT,iSh) = tbl{2,5};
    end
    
    [friedmanClusterSumsSh{iSh}, friedmanMaxAbsClusterSumSh(iSh)] = ...
        rd_clusterSum(friedmanSh(:,iSh), pFriedmanSh(:,iSh)<alpha);
    
    if mod(iSh,10)==0, toc, end
end

% confidence intervals
friedmanCI = prctile(friedmanSh,[2.5 97.5],2);

% cluster sum p-value
friedmanClusterPVal = nnz(friedmanMaxAbsClusterSumSh>friedmanMaxAbsClusterSum)/numel(friedmanMaxAbsClusterSumSh);

% plot
figure
hold on
plot(t, chisqFriedman)
plot(toi, friedmanCI, 'k')
ylabel('Friedman chi-sq')
title('permutation Friedman chi-sq null distribution')

figure
hold on
histogram(friedmanMaxAbsClusterSumSh)
vline(friedmanMaxAbsClusterSum,'k')
ax = gca;
ht = text(friedmanMaxAbsClusterSum*1.1,ax.YLim(2)*0.9,sprintf('p = %.2f', friedmanClusterPVal));
ht.FontSize = 16;
xlabel('Friedman chi-sq max cluster sum')

%% F max cluster window
tIdx = find(t==fvalsMaxClusterWin(1)):find(t==fvalsMaxClusterWin(2));
d = squeeze(mean(data(tIdx,:,:)));

figure
subplot(1,2,1)
bar(d')
xlim([0 nSubjects+1])
legend(condNames)
xlabel('observer')
ylabel('microsaccade rate')
subplot(1,2,2)
errorbar(nanmean(d,2),nanstd(d,0,2)/sqrt(nSubjects - numel(rsidx)),'.','MarkerSize',20)
xlim([0.5 nConds+0.5])
set(gca,'XTick',1:nConds)
set(gca,'XTickLabel',condNames)
rd_supertitle2(sprintf('time window = %d to %d ms', t(tIdx(1)), t(tIdx(end))))

%% permutation t-stat null distribution
nShuffles = 1000;

tvalsSh = []; pTSh = [];
tvalsClusterSumsSh = []; tvalsMaxAbsClusterSumSh = [];
for iSh = 1:nShuffles
    for iSubject = 1:nSubjects
       dataSh(:,:,iSubject) = data(:,randperm(nConds),iSubject);  
    end
    [~, pTSh.nt1(:,iSh), ~, statNT1] = ttest(squeeze(dataSh(:,1,:)-dataSh(:,2,:))');
    [~, pTSh.t2n(:,iSh), ~, statNT2] = ttest(squeeze(dataSh(:,1,:)-dataSh(:,3,:))');
    [~, pTSh.t1t2(:,iSh), ~, statT1T2] = ttest(squeeze(dataSh(:,2,:)-dataSh(:,3,:))');
    
    tvalsSh.t1t2(:,iSh) = statT1T2.tstat;
    tvalsSh.nt1(:,iSh) = statNT1.tstat;
    tvalsSh.t2n(:,iSh) = statNT2.tstat;
    
    if iSh == 1
        compNames = fieldnames(tvalsSh);
        nComps = numel(compNames);
    end
    for iComp = 1:nComps
        comp = compNames{iComp};
        [tvalsClusterSumsSh{iSh,iComp}, tvalsMaxAbsClusterSumSh(iSh,iComp)] = ...
            rd_clusterSum(tvalsSh.(comp)(toiIdx,iSh), pTSh.(comp)(toiIdx,iSh)<alpha);
    end
end

% confidence intervals and cluster sum p-value
compNames = fieldnames(tvalsSh);
nComps = numel(compNames);
for iComp = 1:nComps
    comp = compNames{iComp};
    tvalsCI.(comp) = prctile(tvalsSh.(comp),[2.5 97.5],2); 
    clusterPT.(comp) = nnz(tvalsMaxAbsClusterSumSh(:,iComp)>tvalsMaxAbsClusterSum.(comp))/...
        numel(tvalsMaxAbsClusterSumSh(:,iComp));
end

% plot
figure
for iComp = 1:nComps
    comp = compNames{iComp};
    subplot(nComps,1,iComp)
    hold on
    plot(t, tvals.(comp))
    plot(t, tvalsCI.(comp), 'k')
    title(comp)
%     plotEventTimes([-1000 0 250 750])
end
xlabel('time (ms)')

figure
for iComp = 1:nComps
    comp = compNames{iComp};
    subplot(1,nComps,iComp)
    hold on
    histogram(tvalsMaxAbsClusterSumSh(:,iComp))
    vline(tvalsMaxAbsClusterSum.(comp),'k')
    ax = gca;
    ht = text(tvalsMaxAbsClusterSum.(comp)*1.1,ax.YLim(2)*0.9,sprintf('p = %.3f', clusterPT.(comp)));
    ht.FontSize = 16;
    xlabel('t stat max cluster sum')
    title(comp)
end

%% t max cluster window
tIdx = find(t==tvalsMaxClusterWin.nt1(1)):find(t==tvalsMaxClusterWin.nt1(2));
d = squeeze(mean(data(tIdx,:,:)));

figure
subplot(1,2,1)
bar(d')
xlim([0 nSubjects+1])
legend(condNames)
xlabel('observer')
ylabel('microsaccade rate')
subplot(1,2,2)
errorbar(nanmean(d,2),nanstd(d,0,2)/sqrt(nSubjects-numel(rsidx)),'.','MarkerSize',20)
xlim([0.5 nConds+0.5])
set(gca,'XTick',1:nConds)
set(gca,'XTickLabel',condNames)
rd_supertitle2(sprintf('time window = %d to %d ms', t(tIdx(1)), t(tIdx(end))))

%% bin-based analysis
alpha = 0.05;
binSize = 100;
binStarts = -500:binSize:-binSize; % -499?
binEnds = binStarts+binSize-1;
nBins = numel(binStarts);

% check for normality
hLillieBins = []; hJBBins = [];
for iBin = 1:nBins
    idx = find(t==binStarts(iBin)):find(t==binEnds(iBin));
    vals = squeeze(mean(data(idx,:,:),1))';
    hLillieBins(iBin) = lillietest(vals(:));
    hJBBins(iBin) = jbtest(vals(:));
end

% F test
fvalsBins = []; pFBins = [];
for iBin = 1:nBins
    idx = find(t==binStarts(iBin)):find(t==binEnds(iBin));
    vals = squeeze(mean(data(idx,:,:),1))';
    [fvalsBins(iBin), pFBins(iBin)] = rd_rmANOVA(vals, conds, {'Attention'}, nConds);
end

% Friedman test
friedmanBins = []; pFriedmanBins = [];
for iBin = 1:nBins
    idx = find(t==binStarts(iBin)):find(t==binEnds(iBin));
    vals = squeeze(mean(data(idx,:,:),1))';
    [pFriedmanBins(iBin),tbl] = friedman(vals,1,'off');
    friedmanBins(iBin) = tbl{2,5};
end

% store bin data means
dataBins = [];
for iBin = 1:nBins
    idx = find(t==binStarts(iBin)):find(t==binEnds(iBin));
    dataBins(:,:,iBin) = squeeze(mean(data(idx,:,:),1))';
end
groupMeanBins = squeeze(mean(dataBins,1));
groupSteBins = squeeze(std(dataBins,0,1));

dataBinsNorm = normalizeDC(shiftdim(dataBins,1));
groupMeanBinsNorm = nanmean(dataBinsNorm,3);
groupSteBinsNorm = nanstd(dataBinsNorm,0,3)./sqrt(nSubjects-numel(rsidx));

% for iBin = 1:nBins
%     idx = find(t==binStarts(iBin)):find(t==binEnds(iBin));
%     for iCond = 2:nConds %2
%         if iCond==nConds
%             pairCond = 1; %2
%         else
%             pairCond = 1; %iCond+1; %1
%         end
%         str = sprintf('%s%s', condNames{iCond}, condNames{pairCond});
%         vals = squeeze(mean(data(idx,iCond,:)-data(idx,pairCond,:),1))';
%         [hTBins.(str)(iBin), pTBins.(str)(iBin), ~, stat] = ttest(vals);
%         tvalsBins.(str)(iBin) = stat.tstat;
%     end
% end

% t-test
[tvalsBins, pTBins, hTBins] = deal([]);
switch experiment
    case {'e0','e3','e0e3','e0e3e5','e0_nothresh','e3_nothresh','e0e3e5_nothresh'}
        for iBin = 1:nBins
            idx = find(t==binStarts(iBin)):find(t==binEnds(iBin));
            for iCond = 1:nConds
                if iCond==nConds
                    pairCond = 1;
                else
                    pairCond = iCond+1;
                end
                str = sprintf('%s%s', condNames{iCond}, condNames{pairCond});
                vals = squeeze(mean(data(idx,iCond,:)-data(idx,pairCond,:),1))';
                [hTBins.(str)(iBin), pTBins.(str)(iBin), ~, stat] = ttest(vals);
                tvalsBins.(str)(iBin) = stat.tstat;
            end
        end
    case {'e5','e5_nothresh'}
        for iBin = 1:nBins
            idx = find(t==binStarts(iBin)):find(t==binEnds(iBin));
            for iCond = 2:nConds
                if iCond==nConds
                    pairCond = 2;
                else
                    pairCond = iCond+1;
                end
                str = sprintf('%s%s', condNames{iCond}, condNames{pairCond});
                vals = squeeze(mean(data(idx,iCond,:)-data(idx,pairCond,:),1))';
                [hTBins.(str)(iBin), pTBins.(str)(iBin), ~, stat] = ttest(vals);
                tvalsBins.(str)(iBin) = stat.tstat;
            end
            for iCond = 2:nConds
                pairCond = 1;
                str = sprintf('%s%s', condNames{iCond}, condNames{pairCond});
                vals = squeeze(mean(data(idx,iCond,:)-data(idx,pairCond,:),1))';
                [hTBins.(str)(iBin), pTBins.(str)(iBin), ~, stat] = ttest(vals);
                tvalsBins.(str)(iBin) = stat.tstat;
            end
        end
end

% plots
m = fields(tvalsBins);
nM = numel(m);

figure
subplot(2,1,1)
hold on
plot(binStarts, pFBins, '.', 'MarkerSize', 20)
plot([binStarts(1) binEnds(end)], [alpha alpha], '--k')
legend('F')
subplot(2,1,2)
hold on
for iM = 1:nM
    plot(binStarts, pTBins.(m{iM}), '.', 'MarkerSize', 20)
end
plot([binStarts(1) binEnds(end)], [alpha alpha], '--k')
legend(m)
xlabel('time (ms)')
ylabel('p-value')

figure
errorbar(repmat(binStarts,nConds,1)', groupMeanBins', groupSteBins','LineStyle','none','LineWidth',2)
xlim(binStarts([1 end]) + [-binSize/2 binSize/2])
legend(condNames)
xlabel('Bin (ms)')
ylabel('Microsaccade rate (Hz)')

figure
errorbar(repmat(binStarts,nConds,1)', groupMeanBinsNorm', groupSteBinsNorm','LineStyle','none','LineWidth',2)
xlim(binStarts([1 end]) + [-binSize/2 binSize/2])
legend(condNames)
xlabel('Bin (ms)')
ylabel('Microsaccade rate (Hz)')

%% save
if saveAnalysis
    save(saveFileName, 'subjects', 'condNames', 't', 'eventTimes', 'nSamples', 'toi',...
        'data', 'dataMean', 'dataSte', 'clusters','clusterData','clusterMean','clusterSte')
end

