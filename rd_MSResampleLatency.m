% rd_resamplePrebound.m

load('data/e5/bonneh_post.mat')

trialsData = lastOnset;
trialsMean = mean(boneh_mat,1);

nConds = numel(condNames);
nTrials = [128 256 256 256];
nSamples = 1000;

m = [];
for iCond = 1:nConds
    m(:,iCond,:) = bootstrp(nSamples, @nanmean, trialsData(1:nTrials(iCond),iCond,:));
end
groupMean = mean(m,3);
ci = prctile(groupMean,[2.5 97.5]);
lb = trialsMean-ci(1,:);
ub = ci(2,:)-trialsMean;

% plot 95% confidence intervals
figure('Position',[125 280 940 425])
ax = subplot(1,2,1);
hold on
for iCond = 1:nConds
    histogram(groupMean(:,iCond))
end
subplot(1,2,2)
hold on
bar(1:nConds, trialsMean)
errbar(1:nConds, trialsMean, lb, ub, 'k')  
set(gca,'XTick',1:nConds)
ylim(ax.XLim)
