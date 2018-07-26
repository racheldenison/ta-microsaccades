% rd_MSWriteDataFile.m

% load data using rd_MSStats.m or rd_MSRebound.m

e0S = load('data/e0/subjects.mat');
e3S = load('data/e3/subjects.mat');
e5S = load('data/e5/subjects.mat');
subjects = [e0S.subjects e3S.subjects e5S.subjects];
exp = [repmat(cellstr('e0'),1,numel(e0S.subjects)) ...
    repmat(cellstr('e3'),1,numel(e3S.subjects)) ...
    repmat(cellstr('e5'),1,numel(e5S.subjects))];
        
expNames = {'e0','e3','e5'};
expNums = [ones(1,numel(e0S.subjects)) ...
    2*ones(1,numel(e3S.subjects)) ...
    3*ones(1,numel(e5S.subjects))];

% expNums = 3*ones(1,numel(subjects));

%% time points
% fileID = fopen('E0E3E5_msrate_-500to0ms_N30.txt','w');
fileID = fopen('E0E3E5_msrate_-750to750ms_N30.txt','w');
fprintf(fileID,'%s %s %s %s %s %s\n','experiment','subject','time','neutral','t1','t2');

for iSubject = 1:nSubjects
    exp = expNames{expNums(iSubject)};
    subject = subjects{iSubject};
%     for iT = 1:numel(t)
%         fprintf(fileID,'%s %s %d %1.4f %1.4f %1.4f\n', exp, subject, t(iT), data(iT,:,iSubject));
%     end
    for iT = 1:nT
        fprintf(fileID,'%s %s %d %1.4f %1.4f %1.4f\n', exp, subject, t(toiIdx(iT)), data(toiIdx(iT),:,iSubject));
    end
end
    
fclose(fileID);

%% time bins
fileID = fopen('E0E3E5_msrate_-500to0ms_bin500ms_N30.txt','w');
fprintf(fileID,'%s %s %s %s %s %s %s\n','experiment','subject','bin_start','bin_end','neutral','t1','t2');

for iSubject = 1:nSubjects
    exp = expNames{expNums(iSubject)};
    subject = subjects{iSubject};
    for iBin = 1:size(dataBins,3)
        fprintf(fileID,'%s %s %d %d %1.4f %1.4f %1.4f\n', ...
            exp, subject, binStarts(iBin), binEnds(iBin), dataBins(iSubject,:,iBin));
    end
end
    
fclose(fileID);

%% rebound
% fileID = fopen('E0E3E5_msRebound2Latency_N30_dataZ.txt','w');
fileID = fopen('E5_msRebound2Latency_N9_dataZ.txt','w');
% fprintf(fileID,'%s %s %s %s %s\n','experiment','subject','neutral','t1','t2');
fprintf(fileID,'%s %s %s %s %s %s\n','experiment','subject','neutral','t1','t2','t3');

for iSubject = 1:nSubjects
    exp = expNames{expNums(iSubject)};
    subject = subjects{iSubject};
%     fprintf(fileID,'%s %s %1.4f %1.4f %1.4f\n', exp, subject, dataZ(iSubject,:));
    fprintf(fileID,'%s %s %1.4f %1.4f %1.4f %1.4f\n', exp, subject, dataZ(iSubject,:));
end
    
fclose(fileID);

%% clusters
fileID = fopen('E0E3E5_msrate_-500to0ms_clusterBetaCueNeutral_N30.txt','w');
fprintf(fileID,'%s %s %s %s %s %s\n','cluster','experiment','subject','neutral','t1','t2');

for iCluster = 1:size(clusterData,2)
    for iSubject = 1:nSubjects
        exp = expNames{expNums(iSubject)};
        subject = subjects{iSubject};
        fprintf(fileID,'%d %s %s %1.4f %1.4f %1.4f\n', iCluster, exp, subject, clusterData(:,iCluster,iSubject));
    end
end
    
fclose(fileID);

%% rebound vs. behavior
fileID = fopen('E0E3E5_msRebound2LatencyVBehav_N30.txt','w');
% fileID = fopen('E5_msPrebound2Latency_N9.txt','w');
fprintf(fileID,'%s %s %s %s %s %s %s\n','experiment','subject','target','accuracy','neutral','t1','t2');
% fprintf(fileID,'%s %s %s %s %s %s %s\n','experiment','subject','target','accuracy','neutral','t1','t2','t3');

data = groupDataT;
targetNames = condNames(2:end);

for iT = 1:nTargets
    targetName = targetNames{iT};
    for iAcc = 1:nAccs
        accName = accNames{iAcc};
        for iSubject = 1:nSubjects
            exp = expNames{expNums(iSubject)};
            subject = subjects{iSubject};
            fprintf(fileID,'%s %s %s %s %1.4f %1.4f %1.4f\n', ...
                exp, subject, targetName, accName, data(:,iAcc,iT,iSubject));
            % fprintf(fileID,'%s %s %d %1.4f %1.4f %1.4f %1.4f\n', ...
            %   exp, subject, targetName, accName, data(iSubject,:));
        end
    end
end
    
fclose(fileID);

%% rebound vs. behavior - latency bins
fileID = fopen('E5_msRebound2LatencyVBehavBin_bin200ms_N9_c2500ms.txt','w');
fprintf(fileID,'%s %s %s %s %s %s\n','experiment','subject','target','bin_start','bin_or_mean','acc');

data = accBin;
data0 = accTBin;

for iT = 1:nTargets
    targetName = targetNames{iT};
    for iBin = 1:nBins
        binStart = binStarts(iBin);
        for iSubject = 1:nSubjects
            expName = exp{iSubject}; % exp is from rd_MSReboundVBehavBin.m
            subject = subjects{iSubject};
            fprintf(fileID,'%s %s %s %d bin %1.4f\n', ...
                expName, subject, targetName, binStart, data(iBin,iT,iSubject));
            fprintf(fileID,'%s %s %s %d mean %1.4f\n', ...
                expName, subject, targetName, binStart, data0(iBin,iT,iSubject));
        end
    end
end
    
fclose(fileID);

%% rebound vs. behavior - latency bins - difference from mean
fileID = fopen('E0E3E5_msPrebound2LatencyVBehavBinDiff_t1t2_bin100ms_step10ms_N30.txt','w');
fprintf(fileID,'%s %s %s %s %s\n','experiment','subject','target','bin_start','acc_diff');

data = accBinDiff;

for iT = 1:nTargets
    targetName = targetNames{iT};
    for iBin = 1:nBins
        binStart = binStarts(iBin);
        for iSubject = 1:nSubjects
            expName = exp{iSubject}; % exp is from rd_MSReboundVBehavBin.m
            subject = subjects{iSubject};
            fprintf(fileID,'%s %s %s %d %1.4f\n', ...
                expName, subject, targetName, binStart, data(iBin,iT,iSubject));
        end
    end
end
    
fclose(fileID);

%% ms near target acc
fileID = fopen('E5_msNearTarget_accNT_-100to0_N9.txt','w');
fprintf(fileID,'%s %s %s %s %s\n','experiment','subject','target','ms_in_window','acc');

data = accNT;

for iT = 1:nTargets
    targetName = targetNames{iT};
    for iH = 1:numel(hs)
        h = hs(iH);
        for iSubject = 1:nSubjects
            expName = exp{iSubject}; 
            subject = subjects{iSubject};
            fprintf(fileID,'%s %s %s %d %1.4f\n', ...
                expName, subject, targetName, h, data(iH,iT,iSubject));
        end
    end
end

%% number of ms near target
fileID = fopen('E0E3E5_msNearTarget_msNT_-200to0_N30.txt','w');
fprintf(fileID,'%s %s %s %s %s\n','experiment','subject','target','cue','ms_prop');

data = msNT;
cueTypeNames = {'neutral','target'};

for iT = 1:nTargets
    targetName = targetNames{iT};
    for iC = 1:2
        cueType = cueTypeNames{iC};
        for iSubject = 1:nSubjects
            expName = exp{iSubject}; 
            subject = subjects{iSubject};
            fprintf(fileID,'%s %s %s %s %1.4f\n', ...
                expName, subject, targetName, cueType, data(iC,iT,iSubject));
        end
    end
end

%% accuracy for MS subjects
fileID = fopen('E0E3E5_msBehavAccuracyNormalized_N30.txt','w');
fprintf(fileID,'%s %s %s %s %s\n','experiment','subject','target','validity','acc');

validityNames = {'valid','neutral','invalid'};
nV = numel(validityNames);
nTargets = 2;
nSubjects = numel(subjectInits);

for iT = 1:nTargets
    targetName = targetNames{iT};
    data = b.accData{iT};
    for iV = 1:nV
        validity = validityNames{iV};
        for iSubject = 1:nSubjects
            expName = exp{iSubject}; % exp is from rd_MSReboundVBehavBin.m
            subject = subjectInits{iSubject}(1:2);
            fprintf(fileID,'%s %s %s %s %1.4f\n', ...
                expName, subject, targetName, validity, data(iV,iSubject));
        end
    end
end
    
fclose(fileID);
