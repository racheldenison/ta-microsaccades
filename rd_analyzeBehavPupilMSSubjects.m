function [results, subjectInits, screenInfo] = rd_analyzeBehavPupilMSSubjects(expt, modality, plotFigs)

%% setup
if nargin==0
    expt = 'E0E3E5'; %'E0','E3','E0E3'
    modality = 'MS'; %'pupil','MS'
    plotFigs = 1;
end

switch modality
    case 'pupil'
        subjectInitsE0 = {'ma','ad','bl','ec','ty','zw','hl','rd','jp'}; % Pupil
        subjectInitsE3 = {'blE3','caE3','ecE3','enE3','ewE3','idE3','jlE3','jxE3','ldE3','mlE3','rdE3','sjE3'}; % Pupil
    case 'MS'
        subjectInitsE0 = {'ma' 'ad' 'bl' 'ec' 'ty' 'vp' 'zw' 'hl' 'rd' 'jp'}; % MS
        subjectInitsE3 = {'blE3' 'caE3' 'ecE3' 'enE3' 'ewE3' 'jlE3' 'jxE3' 'ldE3' 'mlE3' 'rdE3' 'sjE3'}; % MS
        subjectInitsE5 = {'dsE5' 'gbE5' 'htE5' 'ikE5' 'jgE5' 'jpE5' 'rdE5' 'xwE5' 'yzE5'}; % MS
    otherwise
        error('modality not recognized')
end         

switch expt
    case 'E0'
        subjectInits = subjectInitsE0;
    case 'E3'
        subjectInits = subjectInitsE3;
    case 'E5'
        subjectInits = subjectInitsE5;
    case 'E0E3'
        subjectInits = [subjectInitsE0 subjectInitsE3];
    case 'E0E3E5'
        subjectInits = [subjectInitsE0 subjectInitsE3 subjectInitsE5];
    otherwise
        error('expt not recognized')
end
tilt = '*';
contrast = '*'; % '64'; % 
contrastIdx = 1; % only plot one contrast at a time
soa1 = 1000;
soa2 = 1250;

run = 9;

normalizeNeutral = 1;
normalizeMorey = 1;

saveFigs = 0;

nSubjects = numel(subjectInits);

%% Get data
for iSubject = 1:nSubjects 
    subjectInit = subjectInits{iSubject};
    
    if strfind(subjectInit, 'E3')
        exptName = 'a1';
        exptDir = 'E3_adjust';
    elseif strfind(subjectInit, 'E5')
        exptName = 'cbD15';
        exptDir = 'E5_T3_cbD15';
    else
        switch subjectInit
            case 'jp'
                exptName = 'cbD15';
                exptDir = 'pilot';
            case {'hl','rd'}
                exptName = 'cbD6';
                exptDir = 'E2_SOA_cbD6';
            otherwise
                exptName = 'cb';
                exptDir = 'E0_cb';
        end
    end
    dataDir = sprintf('%s/%s/%s', pathToExpt('data'), exptDir, subjectInit(1:2));
    
    exptStr = sprintf('%s_%s_soa%d-%d*', exptName, contrast, soa1, soa2);
    subjectID = sprintf('%s*_%s', subjectInit(1:2), exptStr);
    
    % load data from a given soa
    if strfind(subjectInit, 'E5')
        dataFile = dir(sprintf('%s/%s_%s_*_run01WW_TemporalAttention3Targets*', dataDir, subjectInit(1:2), exptName));
    elseif strcmp(subjectInit,'jp')
        dataFile = dir(sprintf('%s/%s_none_run%02d*', dataDir, subjectID, run));
    else
        dataFile = dir(sprintf('%s/%s_run%02d_T*', dataDir, subjectID, run));
    end
    if numel(dataFile)~=1
        sprintf('%s/%s_run%02d*', dataDir, subjectID, run)
        error('more or fewer than one matching data file')
    else
        load(sprintf('%s/%s', dataDir, dataFile.name))
        
        % store a p from E0
        if strcmp(exptDir, 'E0_cb')
            pE0 = expt.p;
        end
        
        % store screen info
        screenInfo.size(iSubject,:) = expt.p.screenSize;
        screenInfo.res(iSubject,:) = expt.p.screenRes;
        screenInfo.viewDist(iSubject,1) = expt.p.viewDist;
        screenInfo.eyeSlack(iSubject,1) = expt.p.eyeSlack;
    end
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% if you want to reanalyze, do it here %%%
%     T1T2Axis = 'same';
%     [expt results] = rd_analyzeTemporalAttention(expt, 0, 0, 0, 0, T1T2Axis, 0);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % read out the accuracy and rt
    for iEL = 1:2 % early/late
        accData{iEL}(:,iSubject) = results.accMean{iEL}(:,contrastIdx);
        rtData{iEL}(:,iSubject) = results.rtMean{iEL}(:,contrastIdx);
    end
    
    % also gather all the data for all contrasts
    for iContrast = 1:numel(expt.p.targetContrasts)
        for iEL = 1:2 % early/late
            accDataAll{iEL}(:,iSubject,iContrast) = results.accMean{iEL}(:,iContrast);
            rtDataAll{iEL}(:,iSubject,iContrast) = results.rtMean{iEL}(:,iContrast);
            
            % average across contrasts
            accDataC{iEL}(:,iSubject) = mean(results.accMean{iEL},2);
            rtDataC{iEL}(:,iSubject) = mean(results.rtMean{iEL},2);
        end
    end
end

%% Normalize and summarize
if normalizeNeutral
    neutralValsAcc = (accDataC{1}(3,:) + accDataC{2}(3,:))/2;
    neutralValsRT = (rtDataC{1}(3,:) + rtDataC{2}(3,:))/2;
    for iEL = 1:2
        accDataC{iEL} = accDataC{iEL}./repmat(neutralValsAcc,3,1);
        rtDataC{iEL} = rtDataC{iEL}./repmat(neutralValsRT,3,1);
    end
end

if normalizeMorey
    % use the Morey 2008 correction
    a = cat(3, accDataC{1}, accDataC{2});
    b = normalizeDC(shiftdim(a,2));
    c = shiftdim(b,1);
    accDataC{1} = c(:,:,1);
    accDataC{2} = c(:,:,2);
    
    a = cat(3, rtDataC{1}, rtDataC{2});
    b = normalizeDC(shiftdim(a,2));
    c = shiftdim(b,1);
    rtDataC{1} = c(:,:,1);
    rtDataC{2} = c(:,:,2);
    
    [fixed1 fixed2 N] = size(b);
    M = fixed1*fixed2;
    morey = M/(M-1);
    
    for iEL = 1:2
        accMeanC(:,iEL) = mean(accDataC{iEL},2);
        accSteC(:,iEL) = sqrt(morey*var(accDataC{iEL},0,2)./(nSubjects));
        rtMeanC(:,iEL) = mean(rtDataC{iEL},2);
        rtSteC(:,iEL) = sqrt(morey*var(rtDataC{iEL},0,2)./(nSubjects));
    end
else
    for iEL = 1:2 % early/late
        accMeanC(:,iEL) = mean(accDataC{iEL},2);
        accSteC(:,iEL) = std(accDataC{iEL},0,2)./sqrt(nSubjects);
        rtMeanC(:,iEL) = mean(rtDataC{iEL},2);
        rtSteC(:,iEL) = std(rtDataC{iEL},0,2)./sqrt(nSubjects);
    end
end

%% Store results
cvOrder = [1 3 2]; % [valid neutral invalid]
for iEL = 1:2
    results.accData{iEL} = accDataC{iEL}(cvOrder,:);
    results.rtData{iEL} = rtDataC{iEL}(cvOrder,:);
end
results.accMean = accMeanC(cvOrder,:);
results.accSte = accSteC(cvOrder,:);
results.rtMean = rtMeanC(cvOrder,:);
results.rtSte = rtSteC(cvOrder,:);

%% Let p come from E0 if possible
if exist('pE0','var')
    p = pE0;
else
    p = expt.p;
end

%% Plot
if plotFigs
    nCV = numel(p.cueValidity);
    intervalNames = {'T1','T2'};
    cueNames = {'valid','invalid','neutral'};
    if normalizeNeutral
        accLims = [0 2];
        rtLims = [0 2];
    else
        accLims = [0.5 0.9];
        rtLims = [0.9 1.6];
    end
    xlims = [0 nSubjects+1];
    axTitle = '';
    [y, idx] = sort(p.cueValidity,2,'descend');
    
    % indiv subjects
    figure
    for iRI = 1:numel(p.respInterval)
        subplot(1,numel(p.respInterval),iRI);
        hold on
        
        p1 = bar(repmat((1:nSubjects)',1,nCV),...
            accDataC{iRI}(idx,:)');
        
        set(gca,'XTick',1:nSubjects)
        xlabel('subject')
        ylabel('acc')
        legend(p1, cueNames(idx),'location','best')
        title(intervalNames{iRI})
        xlim(xlims)
        ylim(accLims)
        rd_supertitle(sprintf('N=%d', nSubjects));
        rd_raiseAxis(gca);
    end
    
    figure
    for iRI = 1:numel(p.respInterval)
        subplot(1,numel(p.respInterval),iRI);
        
        p1 = bar(repmat((1:nSubjects)',1,nCV),...
            rtDataC{iRI}(idx,:)');
        
        set(gca,'XTick',1:nSubjects)
        xlabel('subject')
        ylabel('rt')
        legend(cueNames(idx),'location','best')
        title(intervalNames{iRI})
        xlim(xlims)
        ylim(rtLims)
        box off
        rd_supertitle(sprintf('N=%d', nSubjects));
        rd_raiseAxis(gca);
    end
    
    % group
    figure
    for iRI = 1:numel(p.respInterval)
        subplot(1,numel(p.respInterval),iRI);
        hold on
        
        b1 = bar(1:nCV, accMeanC(idx,iRI),'FaceColor',[.5 .5 .5]);
        p1 = errorbar(1:nCV, accMeanC(idx,iRI)', accSteC(idx,iRI)','k','LineStyle','none');
        
        set(gca,'XTick',1:nCV)
        set(gca,'XTickLabel', cueNames(idx))
        xlabel('cue validity')
        ylabel('acc')
        title(intervalNames{iRI})
        ylim(accLims)
        box off
        rd_supertitle(sprintf('N=%d', nSubjects));
        rd_raiseAxis(gca);
    end
    
    figure
    for iRI = 1:numel(p.respInterval)
        subplot(1,numel(p.respInterval),iRI);
        hold on
        
        b1 = bar(1:nCV, rtMeanC(idx,iRI),'FaceColor',[.5 .5 .5]);
        p1 = errorbar(1:nCV, rtMeanC(idx,iRI)', rtSteC(idx,iRI)','k','LineStyle','none');
        
        set(gca,'XTick',1:nCV)
        set(gca,'XTickLabel', cueNames(idx))
        xlabel('cue validity')
        ylabel('rt')
        title(intervalNames{iRI})
        ylim(rtLims);
        box off
        rd_supertitle(sprintf('N=%d', nSubjects));
        rd_raiseAxis(gca);
    end
end

%% Reshape data for output to R
acc1 = reshape(accDataC{1}', nSubjects*3, 1);
acc2 = reshape(accDataC{2}', nSubjects*3, 1);
rt1 = reshape(rtDataC{1}', nSubjects*3, 1);
rt2 = reshape(rtDataC{2}', nSubjects*3, 1);

acc_all = [acc1; acc2];
rt_all = [rt1; rt2];

%% Quick stats
for iEL = 1:2
    fprintf('T%d\n',iEL)
    vals = rtDataC{iEL};
    [hvi pvi cvi svi] = ttest(vals(1,:),vals(2,:));
    [hvn pvn cvn svn] = ttest(vals(1,:),vals(3,:));
    [hni pni cni sni] = ttest(vals(2,:),vals(3,:));
    fprintf('valid vs. invalid, t(%d) = %1.3f, p = %1.4f\n', svi.df, svi.tstat, pvi)
    fprintf('valid vs. neutral, t(%d) = %1.3f, p = %1.4f\n', svn.df, svn.tstat, pvn)
    fprintf('neutral vs. invalid, t(%d) = %1.3f, p = %1.4f\n\n', sni.df, sni.tstat, pni)
end

%% Effect size
% calculate observed pairwise differences
% ACC
% collapsing across contrast
for iT = 1:2
    accDataCP(1,:,iT) = accDataC{iT}(1,:) - accDataC{iT}(2,:); % VI
    accDataCP(2,:,iT) = accDataC{iT}(1,:) - accDataC{iT}(3,:); % VN
    accDataCP(3,:,iT) = accDataC{iT}(3,:) - accDataC{iT}(2,:); % NI
end

accdCP = mean(accDataCP,2)./std(accDataCP,0,2);

% % collapsing across target and contrast
% accDataCTP(1,:) = accDataCT(1,:) - accDataCT(2,:); % VI
% accDataCTP(2,:) = accDataCT(1,:) - accDataCT(3,:); % VN
% accDataCTP(3,:) = accDataCT(3,:) - accDataCT(2,:); % NI

% accdCTP = mean(accDataCTP,2)./std(accDataCTP,0,2);

% RT
% collapsing across contrast
for iT = 1:2
    rtDataCP(1,:,iT) = rtDataC{iT}(1,:) - rtDataC{iT}(2,:); % VI
    rtDataCP(2,:,iT) = rtDataC{iT}(1,:) - rtDataC{iT}(3,:); % VN
    rtDataCP(3,:,iT) = rtDataC{iT}(3,:) - rtDataC{iT}(2,:); % NI
end

rtdCP = mean(rtDataCP,2)./std(rtDataCP,0,2);

% % collapsing across target and contrast
% rtDataCTP(1,:) = rtDataCT(1,:) - rtDataCT(2,:); % VI
% rtDataCTP(2,:) = rtDataCT(1,:) - rtDataCT(3,:); % VN
% rtDataCTP(3,:) = rtDataCT(3,:) - rtDataCT(2,:); % NI
% 
% rtdCTP = mean(rtDataCTP,2)./std(rtDataCTP,0,2);

% R: pwr.t.test(d = 1.2335, sig.level = .05, power = .8, type = "paired")
% http://www.statmethods.net/stats/power.html
% http://powerandsamplesize.com/Calculators/Test-1-Mean/1-Sample-Equality

