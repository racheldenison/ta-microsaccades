% rd_MSResampleKSTest.m

% load trialsDataZScore from rd_MSRebound.m

% real data
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

% shuffle conditions
nShuffles = 10000;
for iSh = 1:nShuffles
    for iSubject = 1:nSubjects
        condIdx = randperm(nConds);
        trialsDataZScoreSh(:,:,iSubject,iSh) = trialsDataZScore(:,condIdx,iSubject); 
        
        for iCond = 1:nConds
            vals = trialsDataZScoreSh(:,iCond,iSubject,iSh);
            valsn = vals(~isnan(vals));
            idx = randi(numel(valsn),n,1); % bootstrap sample
            trialsDataZScoreSampledSh(:,iCond,iSubject,iSh) = valsn(idx);
        end
    end
end

%% KS test
switch experiment
    case {'e0','e3','e0e3e5','e0_nothresh','e3_nothresh','e0e3e5_nothresh'}
        for iSh = 1:nShuffles
            for iCond = 1:nConds
                if iCond==nConds
                    pairCond = 1;
                else
                    pairCond = iCond+1;
                end
                str = sprintf('%s%s', condNames{iCond}, condNames{pairCond});
                
                vals1 = squeeze(trialsDataZScoreSampledSh(:,iCond,:,iSh));
                vals2 = squeeze(trialsDataZScoreSampledSh(:,pairCond,:,iSh));
                [hKSSh.(str)(iSh), pKSSh.(str)(iSh), statKSSh.(str)(iSh)] = kstest2(vals1(:), vals2(:));
            end
        end
    case {'e5','e5_nothresh'}
        for iSh = 1:nShuffles
            for iCond = 2:nConds
                if iCond==nConds
                    pairCond = 2;
                else
                    pairCond = iCond+1;
                end
                str = sprintf('%s%s', condNames{iCond}, condNames{pairCond});
                
                vals1 = squeeze(trialsDataZScoreSampledSh(:,iCond,:,iSh));
                vals2 = squeeze(trialsDataZScoreSampledSh(:,pairCond,:,iSh));
                [hKSSh.(str)(iSh), pKSSh.(str)(iSh), statKSSh.(str)(iSh)] = kstest2(vals1(:), vals2(:));
            end
            for iCond = 2:nConds
                pairCond = 1;
                str = sprintf('%s%s', condNames{iCond}, condNames{pairCond});
                
                vals1 = squeeze(trialsDataZScoreSampledSh(:,iCond,:,iSh));
                vals2 = squeeze(trialsDataZScoreSampledSh(:,pairCond,:,iSh));
                [hKSSh.(str)(iSh), pKSSh.(str)(iSh), statKSSh.(str)(iSh)] = kstest2(vals1(:), vals2(:));
            end
        end
end

%% calcualte thresholds
alpha = 0.05;
for iCond = 1:nConds
    if iCond==nConds
        pairCond = 1;
    else
        pairCond = iCond+1;
    end
    str = sprintf('%s%s', condNames{iCond}, condNames{pairCond});
    
    thresh95.(str) = prctile(statKSSh.(str),(1-alpha)*100);
    permPKS.(str) = mean(statKSSh.(str) > statKS.(str));
end

