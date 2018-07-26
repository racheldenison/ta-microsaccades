% rd_MSRegression.m

% see also glmfit

experiment = 'e0';

load(sprintf('data/%s/onsets.mat',experiment))
load(sprintf('data/%s/behav_data.mat',experiment))

nSubjects = numel(onsetsAll);
t = 1:size(onsetsAll{1},2);

nConds = numel(condNames);

switch experiment
    case {'e0','e3'}
        nTrialsPerCueCond = [128 256 256];
        targetNames = {'t1','t2'};
    otherwise
        error('experiment not recognized')
end

nTargets = numel(targetNames);

%% reorganize accData and targetCond
for iS = 1:nSubjects
    accDataAll{iS} = [];
    targetCondAll{iS} = [];
    for iCond = 1:nConds
        accDataAll{iS} = [accDataAll{iS}; accData(1:nTrialsPerCueCond(iCond),iCond,iS)];
        targetCondAll{iS} = [targetCondAll{iS}; respCueCond(1:nTrialsPerCueCond(iCond),iCond,iS)];
    end
end

%% regression on accuracy
B = [];
for iS = 1:nSubjects
    targetCond = targetCondAll{iS};
    for iT = 1:nTargets
        onsets = onsetsAll{iS}(targetCond==iT,:);
        acc = accDataAll{iS}(targetCond==iT,:);
        
        % downsample
        nDS = 125; %125
        onsetsDS = zeros(size(onsets(:,1:nDS:end)));
        for iD = 1:nDS
            onsetsDS = onsetsDS + onsets(:,iD:nDS:end);
        end
        
%         B(:,iT,iS) = glmfit(onsetsDS,acc,'binomial');
        B(:,iT,iS) = mnrfit(onsetsDS,categorical(acc));
    end
end

BMean = nanmean(B,3);

%% logistic regression on cue condition
B = [];
for iS = 1:nSubjects
    onsets = onsetsAll{iS};
    cond = condsAll{iS};
%     cond = -cond;
    
    % downsample
    nDS = 125;
    onsetsDS = zeros(size(onsets(:,1:nDS:end)));
    for iD = 1:nDS
        onsetsDS = onsetsDS + onsets(:,iD:nDS:end);
    end
    
    B(:,:,iS) = mnrfit(onsetsDS,categorical(cond));
end

BMean = nanmean(B,3);

%% plot
names = targetNames; % 'targetNames' for regression on accuracy

figure
plot(t(1:nDS:end),BMean(2:end,:))
xlabel('time (ms)')
ylabel('regression weight')
legend(names)
box off

figure
hold on
for iS = 1:nSubjects
    plot(t(1:nDS:end),B(2:end,:,iS))
    set(gca,'ColorOrderIndex',1) % reset color order
end
xlabel('time (ms)')
ylabel('regression weight')
legend(names)