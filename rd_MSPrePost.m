% rd_MSPrePost.m

%% settings
analysisName1 = 'bonneh_pre2'; %'bonneh','bonneh_pre','bonneh_pre2','bonneh_post','bonneh_post2'
analysisName2 = 'bonneh_post2';
threshExt = ''; % '', '_nothresh'

expName = sprintf('e0e3e5%s', threshExt);

fileName1L = sprintf('data/%s/LatencyAnalysis_%s.mat', expName, analysisName1);
fileName2L = sprintf('data/%s/LatencyAnalysis_%s.mat', expName, analysisName2);

fileName1R = sprintf('data/%s/RateAnalysis.mat', expName);

latencyEventTime = 1500;

%% load data
E1L = load(fileName1L);
E2L = load(fileName2L);

ER = load(fileName1R);

%% set a few values
condNames = E1L.condNames;
nConds = numel(condNames);
nSubjects = numel(E1L.subjects);

rate = squeeze(mean(mean(ER.data,1),2));

figure
colors = get(gca,'ColorOrder');
close(gcf)
condColors = [.5 .5 .5; colors(1:3,:)]; 

measure = 'dataZ'; % 'data','dataZ','dataM'

%% get data
vals1 = E1L.(measure);
vals2 = E2L.(measure);

%% plot setup
switch measure
    case 'dataZ'
        xlims = [-1 1];
        ylims = [-1 1];
        xlimsShift = [-1.5 1.5];
        ylimsShift = [-1.5 1.5];
    case {'data','dataM'}
        xlims = [600 1400];
        ylims = [1800 2600];
        xlimsShift = [-250 250];
        ylimsShift = [-250 250];
    otherwise
        xlims = 'auto';
        ylims = 'auto';
        xlimsShift = 'auto';
        ylimsShift = 'auto';
end

%% pre vs. post
xgrid = linspace(xlims(1),xlims(2),20);
figure('Position',[200 30 800 350])
for iCond = 1:nConds
    p = polyfit(vals1(:,iCond),vals2(:,iCond),1);
    yfit = polyval(p,xgrid);
    
    subplot(1,nConds,iCond)
    hold on
    plot(vals1(:,iCond), vals2(:,iCond), '.', 'Color', condColors(iCond,:), 'MarkerSize', 20)
    plot(xgrid, yfit, 'k')
    axis square
    xlim(xlims)
    ylim(ylims)
    xlabel('pre')
    ylabel('post')
    title(condNames{iCond})
end

%% pre vs. post averaged across cue cond
vals1CueMean = mean(vals1,2);
vals2CueMean = mean(vals2,2);

p = polyfit(vals1CueMean,vals2CueMean,1);
yfit = polyval(p,xgrid);
controlLine = [xgrid(1) xgrid(end); yfit(1) yfit(end)];

[rho, pval] = corr(vals1CueMean, vals2CueMean);

figure
hold on
plot(vals1CueMean, vals2CueMean, '.k', 'MarkerSize', 20)
plot(xgrid, yfit, 'k')
axis square
xlim(xlims)
ylim(ylims)
xlabel('pre')
ylabel('post')
title('mean across cue conditions')

%% definitions of 'group'
% postCutoff = 2272;
% group = vals2CueMean<postCutoff;
% group = control<=median(control);
load group
load control

%% pre shift vs. post shift (t1 to n)
nT1Shift1 = vals1(:,strcmp(condNames,'n')) - vals1(:,strcmp(condNames,'t1'));
nT1Shift2 = vals2(:,strcmp(condNames,'n')) - vals2(:,strcmp(condNames,'t1'));

nT1ShiftCorr = corr(nT1Shift1, nT1Shift2);

figure
hold on
plot(nT1Shift1, nT1Shift2, '.', 'Color','k', 'MarkerSize',20)
plot(nT1Shift1(group), nT1Shift2(group), 'o', 'Color','k', 'MarkerFaceColor','w', 'MarkerSize',5)
axis square
xlim(xlimsShift)
ylim(ylimsShift)
xlabel('pre shift (N-T1)')
ylabel('post shift (N-T1)')

% individual subjects
figure
bar([nT1Shift1 nT1Shift2])

% low vs. high "control"
figure
subplot(1,2,1)
hold on
bar([mean(nT1Shift1(group)) mean(nT1Shift1(~group))],'FaceColor','k')
errbar([1 2],[mean(nT1Shift1(group)) mean(nT1Shift1(~group))],...
    [std(nT1Shift1(group))/sqrt(nnz(group)) mean(nT1Shift1(~group))/sqrt(nnz(~group))])
set(gca,'XTick',[1 2],'XTickLabel',{'low','high'})
ylabel('shift (N-T1)')
title('pre')
subplot(1,2,2)
hold on
bar([mean(nT1Shift2(group)) mean(nT1Shift2(~group))],'FaceColor','k')
errbar([1 2],[mean(nT1Shift2(group)) mean(nT1Shift2(~group))],...
    [std(nT1Shift2(group))/sqrt(nnz(group)) mean(nT1Shift2(~group))/sqrt(nnz(~group))])
title('post')
set(gca,'XTick',[1 2],'XTickLabel',{'low','high'})

% calculate "control"
% Q = [];
% for iSubject = 1:nSubjects
%     P = [vals1CueMean(iSubject), vals2CueMean(iSubject)]';
%     Q(iSubject,:) = projectPointOntoLine(P, controlLine(:,1), controlLine(:,2));
% end
% controlLims = [Q(Q(:,2)==min(Q(:,2)),:); Q(Q(:,2)==max(Q(:,2)),:)];
% controlLength = norm(controlLims(1,:)-controlLims(2,:));
% control = [];
% for iSubject = 1:nSubjects
%     control(iSubject,1) = norm(Q(iSubject,:)-controlLims(1,:))/controlLength;
% end

% shift vs. "control"
xlims = [0 1];
% ylims = [-1.5 1.5];
xgrid = linspace(xlims(1),xlims(2),20);

figure
subplot(1,2,1)
hold on
p = polyfit(control,nT1Shift1,1);
yfit = polyval(p,xgrid);
plot(control, nT1Shift1, '.', 'Color','k', 'MarkerSize',20)
plot(xgrid, yfit, 'k')
% ylim(ylims)
xlabel('control')
ylabel('shift (N-T1)')
title('pre')
subplot(1,2,2)
hold on
p = polyfit(control,nT1Shift2,1);
yfit = polyval(p,xgrid);
plot(control, nT1Shift2, '.', 'Color','k', 'MarkerSize',20)
plot(xgrid, yfit, 'k')
% ylim(ylims)
title('post')




