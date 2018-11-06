% rd_MSPlotCombinedFigs.m

%% settings
analysisName = 'bonneh_pre2'; %'bonneh','bonneh_pre','bonneh_pre2','bonneh_post','bonneh_post2'
threshExt = ''; % '', '_nothresh'

expName1 = sprintf('e0e3e5%s', threshExt);
expName2 = sprintf('e5%s', threshExt);

fileName1L = sprintf('data/%s/LatencyAnalysis_%s.mat', expName1, analysisName);
fileName2L = sprintf('data/%s/LatencyAnalysis_%s.mat', expName2, analysisName);

fileName1R = sprintf('data/%s/RateAnalysis.mat', expName1);
fileName2R = sprintf('data/%s/RateAnalysis.mat', expName2);

fileName1LBin = sprintf('data/%s/LatencyVBehavBinAnalysis_%s.mat', expName1, analysisName);
fileName2LBin = sprintf('data/%s/LatencyVBehavBinAnalysis_%s.mat', expName2, analysisName);

fileName1LN = sprintf('data/%s/LatencyNearTargetAnalysis.mat', expName1);
fileName2LN = sprintf('data/%s/LatencyNearTargetAnalysis.mat', expName2);

latencyEventTime = 1500;

%% load data
E1L = load(fileName1L);
E2L = load(fileName2L);

E1R = load(fileName1R);
E2R = load(fileName2R);

E1LBin = load(fileName1LBin);
E2LBin = load(fileName2LBin);

E1LN = load(fileName1LN);
E2LN = load(fileName2LN);

%% set a few values
E1L.nConds = numel(E1L.condNames);
E2L.nConds = numel(E2L.condNames);

E1L.nSubjects = numel(E1L.subjects);
E2L.nSubjects = numel(E2L.subjects);

figure
colors = get(gca,'ColorOrder');
close(gcf)
condColors = [.5 .5 .5; colors(1:3,:)]; 
% condColors = [127 127 127; 121 129 187; 250 92 64; 92 160 40]/255; 

%% time series with shaded significant clusters
clusters = E1R.clusters;
dataMean = E1R.dataMean;
nConds = numel(E1R.condNames);
condNames = E2R.condNames;
t = E1R.t;
eventTimes = [-1000 E2R.eventTimes];
eventNames = {'precue','T1','T2','T3'};
nEvents = numel(eventTimes);

cueColors = {
    [.5 .5 .5] 
    [123 146 202]./255
    [242 123 96]./255
    [249 201 110]./255};

targetColors = {
    [0 .45 .74]
    [.85 .325 .1]
    [.93 .69 .125]};

clusterNames = fields(clusters);
clusterColors = {[.8 .8 .8], [251 197 185]./255};
plotClusters = 1;

if plotClusters
    figPos = [150 800 1000 450];
    xlims = [-500 0];
    ylims = [0 1.8];
    lineWidth = 4;
    legendLocation = 'NorthEast';
else
    figPos = [150 800 1200 450];
    xlims = [-1000 1000];
    ylims = [0 2.5];
    lineWidth = 3;
    legendLocation = 'NorthEastOutside';
end

figure('Position', figPos)
hold on
if plotClusters
    for iComp = 1:numel(clusterNames)
        clusterName = clusterNames{iComp};
        vals = clusters.(clusterName);
        for iC = 1:size(vals,1)
            clusterTimes = vals(iC,2:3);
            rectangle('Position',[clusterTimes(1), ylims(1), diff(clusterTimes), diff(ylims)],...
                'FaceColor', clusterColors{iComp}, 'EdgeColor', 'none')
        end
    end
end
for i = 1:nConds
    plot(t, dataMean(:,i), 'Color', cueColors{i},'LineWidth',lineWidth)
end
plot(t, E2R.dataMean(:,4),'Color', cueColors{4}, 'LineWidth',lineWidth)
for i = 1:nConds
    plot(t, dataMean(:,i), 'Color', cueColors{i},'LineWidth',lineWidth)
end
ylim(ylims)
for i = 1:nEvents
    vline(eventTimes(i),'Color','k','LineStyle',':','LineWidth',2.5)
    text(eventTimes(i) + diff(xlims)*.01, ylims(2)-.1, eventNames{i},'FontSize',24)
end
legend(condNames)
xlim(xlims)
xlabel('Time (ms)')
ylabel('Microsaccade rate (Hz)')

legend('neutral','cue T1','cue T2','cue T3','Location',legendLocation)
legend boxoff
set(gca,'LineWidth',1)
set(gca,'TickDir','out')
set(gca,'FontSize',24)


%% mean latency
groupMeanNorm = [E1L.groupMeanNorm(2:3) E2L.groupMeanNorm(4) E1L.groupMeanNorm(1)];
groupMeanNorm = groupMeanNorm - latencyEventTime;
groupSteNorm = [E1L.groupSteNorm(2:3) E2L.groupSteNorm(4) E1L.groupSteNorm(1)];
% condNames = [E1.condNames(2:3) E2.condNames(4) E1.condNames(1)];
condNames = {'T1','T2','T3','neutral'};

figure
hold on
h = errbar(1:4, groupMeanNorm, groupSteNorm,'color',[.7 .7 .7],'LineWidth',4);
plot(1:4, groupMeanNorm, '.k', 'MarkerSize', 40)
plot(3, groupMeanNorm(3), '.', 'Color',[.5 .5 .5], 'MarkerSize', 40)
xlim([0.5 numel(condNames)+0.5])
% ylim([760 820])
set(gca,'XTick',1:numel(condNames))
set(gca,'XTickLabel',condNames)
xlabel('Cue')
ylabel('Latency (ms)')
title(und2space(analysisName))
h(3).Color = [.9 .9 .9];

%% mean latency - no T3
normalizeOption = 'morey'; % 'none','morey'

switch normalizeOption
    case 'morey'
        M = 3; % 3 bars for neutral, T1, T2
        m = M/(M-1);
    otherwise
        m = 1;
end

groupMeanNorm = [E1L.groupMeanNorm(2:3) E1L.groupMeanNorm(1)];
groupMeanNorm = groupMeanNorm - latencyEventTime;
groupSteNorm = [E1L.groupSteNorm(2:3) E1L.groupSteNorm(1)]*sqrt(m);
condNames = {'T1','T2','neutral'};
cols = cueColors([2 3 1]);
barCols = [198 208 231; 250 196 184; 198 198 198]./255;

figure
hold on
% h = errbar(1:3, groupMeanNorm, groupSteNorm,'color',[.7 .7 .7],'LineWidth',4);
% plot(1:3, groupMeanNorm, '.k', 'MarkerSize', 40)
for iCond = 1:numel(condNames)
    h = errbar(iCond, groupMeanNorm(iCond), groupSteNorm(iCond),'color',barCols(iCond,:),'LineWidth',8);
    plot(iCond, groupMeanNorm(iCond), '.', 'color', cols{iCond}, 'MarkerSize', 80)
end
xlim([0.5 numel(condNames)+0.5])
set(gca,'XTick',1:numel(condNames))
set(gca,'XTickLabel',condNames)
xlabel('Cue')
ylabel('Latency (ms)')
set(gca,'LineWidth',1)
set(gca,'TickDir','out')
set(gca,'FontSize',36)

%% mean latency zscore
groupZMean = [E1L.groupZMean(2:3) E2L.groupZMean(4) E1L.groupZMean(1)];
groupZSte = [E1L.groupZSte(2:3) E2L.groupZSte(4) E1L.groupZSte(1)];
condNames = {'T1','T2','T3','neutral'};
cols = cueColors([2 3 4 1]);
barCols = [198 208 231; 250 196 184; 253 231 192; 198 198 198]./255; % 252 220 164

figure
hold on
for iCond = 1:numel(condNames)
    h = errbar(iCond, groupZMean(iCond), groupZSte(iCond),'color',barCols(iCond,:),'LineWidth',8);
    plot(iCond, groupZMean(iCond), '.', 'color', cols{iCond}, 'MarkerSize', 80)
end
xlim([0.5 numel(condNames)+0.5])
% ylim([-.3 .3])
set(gca,'XTick',1:numel(condNames))
set(gca,'XTickLabel',condNames)
xlabel('Cue')
ylabel('Latency (z-score)')
set(gca,'LineWidth',2)
set(gca,'TickDir','out')
set(gca,'FontSize',36)
% pbaspect([.9 1 1])

%% zscore density
% check xgrids
if E1L.xgridz ~= E2L.xgridz
    error('E1L and E2L xgridz are not equal - check plotting')
end
% calculate mean density, median of each density, and median crossing point
trialsDataZScoreDensityMean = cat(2, mean(E1L.trialsDataZScoreDensity,3), ...
    mean(E2L.trialsDataZScoreDensity(:,4,:),3));
e1ZMedian = squeeze(nanmedian(nanmedian(E1L.trialsDataZScore,3)));
e2ZMedian = squeeze(nanmedian(nanmedian(E2L.trialsDataZScore,3)));
zMedian = [e1ZMedian e2ZMedian(4)];
for iCond = 1:numel(zMedian)
    idx = find(abs(E1L.xgridz-zMedian(iCond)) == min(abs(E1L.xgridz-zMedian(iCond))));
    zMedianYVal(iCond) = trialsDataZScoreDensityMean(idx,iCond);
end
transparency = 0;
figure('Position',[150 50 1000 450]*.8)
hold on
for iCond = 1:4
    plot(E1L.xgridz, trialsDataZScoreDensityMean(:,iCond),'color',cueColors{iCond})
end
% E1 neutral
shadedErrorBar(E1L.xgridz, mean(E1L.trialsDataZScoreDensity(:,1,:),3), ...
    std(E1L.trialsDataZScoreDensity(:,1,:),0,3)./sqrt(E1L.nSubjects), ...
    {'Color',cueColors{1}},transparency);
% E2 t3
shadedErrorBar(E2L.xgridz, mean(E2L.trialsDataZScoreDensity(:,4,:),3), ...
    std(E2L.trialsDataZScoreDensity(:,4,:),0,3)./sqrt(E2L.nSubjects), ...
    {'Color',cueColors{4}},transparency);
% E1 t2, t1
for iCond = [3 2]
    shadedErrorBar(E1L.xgridz, mean(E1L.trialsDataZScoreDensity(:,iCond,:),3), ...
        std(E1L.trialsDataZScoreDensity(:,iCond,:),0,3)./sqrt(E1L.nSubjects), ...
        {'Color',cueColors{iCond}},transparency);
end
% vertical lines at median
for iCond = 1:4
    plot([zMedian(iCond) zMedian(iCond)],[0 zMedianYVal(iCond)],'color',cueColors{iCond})
end
xlim([-3 3])
ylim([0 .7])
legend('neutral','cue T1','cue T2','cue T3')
legend boxoff
xlabel('Latency (z-score)')
ylabel('Probability density')
set(gca,'LineWidth',1)
set(gca,'TickDir','out')
% title(und2space([analysisName threshExt]))
set(gca,'FontSize',24)

%% plot histogram for one subject
% trialsData has no normalization
iSubject = 7; % rd e5
ylims = [0 50];
figure
hold on
for iCond = 1:4
    histogram(E2L.trialsData(:,iCond,iSubject),'BinWidth',50,'FaceColor',cueColors{iCond})
end
plot([500 500],ylims,'--k','LineWidth',2)
plot([1500 1500],ylims,'--k','LineWidth',2)
text(530,ylims(2),'precue','FontSize',20,'VerticalAlign','top')
text(1530,ylims(2),'T1','FontSize',20,'VerticalAlign','top')
xlabel('Latency (ms)')
ylabel('Number of trials')
legend('neutral','cue T1','cue T2','cue T3','Location','best')
legend boxoff
xlim([400 1600])
ylim(ylims)
set(gca,'XTick',[500 1000 1500])
set(gca,'XTickLabel',[-1000 -500 0])
set(gca,'TickDir','out','FontSize',20,'LineWidth',1)

%% plot latency bin vs. behavior
binStarts = E1LBin.binStarts - latencyEventTime;
targetTimes = E2LBin.targetTimes - latencyEventTime;
nTargets = numel(targetTimes);
cueDur = 200;
targetDur = 30;
binBoundaries = [binStarts binStarts(end)+E1LBin.binSize];

% diff from mean
xlims = [binStarts(1) binStarts(end)+E1LBin.binSize];
ylims = [-.26 .26];
j = [-10 0 10];
grayBoxColor = [.9 .9 .9];

figure
hold on
% plotboxes([eventTimes' (eventTimes+[cueDur repmat(targetDur,1,nTargets)])'],ylims, grayBoxColor)
plotboxes([targetTimes' (targetTimes+repmat(targetDur,1,nTargets))'],ylims, grayBoxColor)
for iB = 1:numel(binBoundaries)
    plot([binBoundaries(iB) binBoundaries(iB)],ylims,':k')
end
plot(xlims,[0 0],'k')
for iT = 1:nTargets-1
    errbar(binStarts+E1LBin.binSize/2+j(iT),E1LBin.accBinDiffMean(:,iT),E1LBin.accBinDiffSte(:,iT),...
        'Color',targetColors{iT},'LineWidth',2)
    p1(iT) = plot(binStarts+E1LBin.binSize/2+j(iT),E1LBin.accBinDiffMean(:,iT),...
        'Color',targetColors{iT},'Marker','.','MarkerSize',30,'LineWidth',2);
end
iT = 3;
errbar(binStarts+E2LBin.binSize/2+j(iT),E2LBin.accBinDiffMean(:,iT),E2LBin.accBinDiffSte(:,iT),...
    'Color',targetColors{iT},'LineWidth',2)
p1(iT) = plot(binStarts+E2LBin.binSize/2+j(iT),E2LBin.accBinDiffMean(:,iT),...
    'Color',targetColors{iT},'Marker','.','MarkerSize',30,'LineWidth',2);
xlabel('Latency (ms)')
ylabel('\Delta proportion correct')
xlim(xlims)
ylim(ylims)
legend(p1, {'target T1','target T2','target T3'})
legend boxoff
set(gca,'TickDir','out','LineWidth',1,'FontSize',20)
set(gca,'XTick',binBoundaries)
% set(gca,'XTick',binStarts+E1LBin.binSize/2)
if strfind(analysisName, 'pre')
    pbaspect([1 1 1])
else
    pbaspect([6/5 1 1])
    for i = 1:nTargets
        text(targetTimes(i) + 40, ylims(2), E2LBin.targetNames{i},'FontSize',20,'VerticalAlign','top')
    end
end

%% onsets
% or load example_onsets.mat from VSS poster Figures folder
expDir = 'data/e5'; 
A = load(sprintf('%s/analysis_struct', expDir));
microsaccades_data = A.analysis_struct.microsaccades_data;
iSubject = 1; condName = 'EVENT_CUE_1';
onsets = microsaccades_data{iSubject}.(condName).onsets;

nTrials = 40;
plotPrePostShading = 1;
if plotPrePostShading
    rasterColor = [0 0 0];
else
    rasterColor = cueColors{2};
end

figure('Position', [100 100 1000 400]) 
hold on
for iTrial = 1:nTrials
    onsetTimes = find(onsets(iTrial,:)) - latencyEventTime;
    nOnsets = length(onsetTimes);
    if plotPrePostShading
        lastPreT1Onset = onsetTimes(find(onsetTimes<0,1,'last'));
        firstPostT1Onset = onsetTimes(find(onsetTimes>0,1,'first'));
        plot([lastPreT1Onset lastPreT1Onset], [1 1].*iTrial + [-.2 1.1], ...
            'Color', cueColors{4}, 'LineWidth', 7)
        plot([firstPostT1Onset firstPostT1Onset], [1 1].*iTrial + [-.2 1.1], ...
            'Color', cueColors{3}, 'LineWidth', 7)
    end
    for iOnset = 1:nOnsets
        plot([onsetTimes(iOnset) onsetTimes(iOnset)], [1 1].*iTrial + [0 .9], ...
            'Color', rasterColor, 'LineWidth', 3) % cueColors{2}
    end
end
for i = 1:nEvents
    vline(eventTimes(i),'Color','k','LineStyle',':','LineWidth',2)
    text(eventTimes(i) + 30, nTrials+5, eventNames{i},'FontSize',36)
end
vline(eventTimes(end)+500,'Color','k','LineStyle',':','LineWidth',2)
text(eventTimes(end)+500 + 30, nTrials+5, 'resp cue','FontSize',36)
xlim([-1000 2000])
ylim([0 nTrials+5])
xlabel('Time (ms)')
ylabel('Trial')
set(gca,'TickDir','out','LineWidth',1,'FontSize',36)

%% MS near target
targetNames = {'T1','T2','T3'};
nTargets = numel(targetNames);

% accNT
accNTDiffMean = [E1LN.accNTDiffMean; E2LN.accNTDiffMean(3)];
accNTDiffSte = [E1LN.accNTDiffSte; E2LN.accNTDiffSte(3)];

figure
hold on
plot([.7 nTargets+.3],[0 0],'--k')
errorbar(1:nTargets, accNTDiffMean, accNTDiffSte,'.k','MarkerSize',30,'LineWidth',2)
set(gca,'XTick',1:nTargets,'XTickLabel',targetNames)
xlabel('Target')
ylabel('\Delta proportion correct (MS - no MS)')
xlim([.7 nTargets+.3])
ylim([-.2 .2])
pbaspect([.5 1 1])
set(gca,'TickDir','out','LineWidth',1,'FontSize',20)

% msNT
msNTMean = [E1LN.msNTMean E2LN.msNTMean(:,3)];
msNTSte = [E1LN.msNTSte E2LN.msNTSte(:,3)];

figure
hold on
p1 = errorbar(msNTMean', msNTSte','.','MarkerSize',30,'LineWidth',2);
p1(1).Color = cueColors{1};
p1(2).Color = [0 0 0];
set(gca,'XTick',1:nTargets,'XTickLabel',targetNames)
xlabel('Target')
ylabel('Proportion of trials with MS')
legend({'neutral','cue target'},'box','off')
xlim([.7 nTargets+.3])
ylim([0 .04])
pbaspect([.5 1 1])
box off
jitterx(gca);
set(gca,'TickDir','out','LineWidth',1,'FontSize',20)

%% behavior
addpath('/Local/Users/denison/Google Drive/NYU/Projects/Temporal_Attention/Code/Pupil')
addpath('/Local/Users/denison/Google Drive/NYU/Projects/Temporal_Attention/Code/Expt_Scripts/Behav')

bExpt = 'E0E3E5';
modality = 'MS';
b = rd_analyzeBehavPupilMSSubjects(bExpt,modality,0);
bMeasures = {'acc','rt'};

t1colors = {
    [140 76 146]/255
    [164 112 168]/255
    [185 147 189]/255};

t2colors = {
    [32 168 138]/255
    [97 181 156]/255
    [138 196 177]/255};

accLims = [0 1.2];
rtLims = [0 1.2];

nTs = 2;
nVs = 3;
vNamesShort = {'V','N','I'};

saveFigs = 0;

for iM = 1:numel(bMeasures)
    m = bMeasures{iM};
    mMean = strcat(m,'Mean');
    mSte = strcat(m,'Ste');
    ylims = eval(strcat(m,'Lims'));
    if strcmp(m,'acc')
        s = 1;
        ylab = 'Normalized accuracy';
    else
        s = 1;
        ylab = 'Normalized RT';
    end
    
    figure
    for iT = 1:nTs
        if iT==1
            colors = t1colors;
        else
            colors = t2colors;
        end
        subplot(1,nTs,iT);
        hold on
        
        for iV = 1:nVs
        	bar(iV, b.(mMean)(iV,iT)*s,'FaceColor',colors{iV},'EdgeColor','none');
        end
        errbar(1:nVs, b.(mMean)(:,iT)'*s, b.(mSte)(:,iT)'*s,'k','LineWidth',2);
        
        set(gca,'XTick',1:nVs)
        set(gca,'XTickLabel', vNamesShort)
        ylabel(ylab)
        ylim(ylims*s)
        box off
        set(gca,'TickDir','out','LineWidth',1,'FontSize',16)
        pbaspect([.7 1 1])
    end
    
    if saveFigs
        print_pdf(sprintf('~/Desktop/Behav_%s_%s',m, bExpt))
    end
end



