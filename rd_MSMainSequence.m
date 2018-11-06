% rd_MSMainSequence.m

%% settings
experiment = 'e0e3e5'; %'e0','e3','e5','e0e3','e0e3e5'

saveFileName = sprintf('data/%s/MainSequenceAnalysis.mat', experiment);
saveAnalysis = true;

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
        analysis_struct.microsaccades_data = ...
            [e0.analysis_struct.microsaccades_data, e3.analysis_struct.microsaccades_data];
        e0S = load(sprintf('data/%s/subjects.mat', 'e0')); 
        e3S = load(sprintf('data/%s/subjects.mat', 'e3')); 
        subjects = [e0S.subjects e3S.subjects];
    case 'e0e3e5'
        e0 = load('data/e0/analysis_struct4.mat');
        e3 = load('data/e3/analysis_struct-11_subjects.mat');
        e5 = load('data/e5/analysis_struct.mat');
        analysis_struct.microsaccades_data = ...
            [e0.analysis_struct.microsaccades_data, e3.analysis_struct.microsaccades_data, ...
            e5.analysis_struct.microsaccades_data];
        e0S = load(sprintf('data/%s/subjects.mat', 'e0'));
        e3S = load(sprintf('data/%s/subjects.mat', 'e3'));
        e5S = load(sprintf('data/%s/subjects.mat', 'e5'));
        subjects = [e0S.subjects e3S.subjects e5S.subjects];
    otherwise
        load(sprintf('data/%s/analysis_struct.mat', experiment))
        load(sprintf('data/%s/subjects.mat', experiment))
end

%% setup
measures = {'amplitudes','velocities'};
nM = numel(measures);

dataStruct = analysis_struct.microsaccades_data;
nSubjects = numel(dataStruct);

conds = fields(dataStruct{1});
nConds = numel(conds);

%% organize data
data = [];
for iSubject = 1:nSubjects
    for iCond = 1:nConds
        for iM = 1:nM
            measure = measures{iM};
            data{iCond,iSubject}.(measure) = dataStruct{iSubject}.(conds{iCond}).(measure);
        end
    end
end

%% aggregate data by subject
for iSubject = 1:nSubjects
    for iM = 1:nM
        measure = measures{iM};
        groupData.(measure){iSubject} = [];
        for iCond = 1:nConds
            groupData.(measure){iSubject} = [groupData.(measure){iSubject}, data{iCond,iSubject}.(measure)];
            groupDataLog.(measure){iSubject} = log(groupData.(measure){iSubject});
        end
    end
end

%% calculate average amplitude for each subject
for iSubject = 1:nSubjects
    for iM = 1:nM
        measure = measures{iM};
        groupData.means.(measure)(iSubject) = mean(groupData.amplitudes{iSubject});
    end
end

%% calculate correlation for each subject
ampAll = []; velAll = [];
for iSubject = 1:nSubjects
    ampVelCorr(iSubject) = corr(groupDataLog.amplitudes{iSubject}', groupDataLog.velocities{iSubject}');
    ampAll = [ampAll, groupDataLog.amplitudes{iSubject}];
    velAll = [velAll, groupDataLog.velocities{iSubject}];
end

ampVelCorrMean = mean(ampVelCorr);
ampVelCorrStd = std(ampVelCorr);

ampVelAllCorr = corr(ampAll',velAll');

%% plot main sequence
figure
hold on
for iSubject = 1:nSubjects
    plot(groupDataLog.amplitudes{iSubject}, groupDataLog.velocities{iSubject}, '.')
%     xlim([0 1])
%     ylim([0 160])
%     pause(1)
end
xlabel('log amplitude')
ylabel('log velocity')

%% save
if saveAnalysis
    save(saveFileName, 'experiment','subjects','measures','conds',...
        'groupData','groupDataLog','ampVelCorr','ampAll','velAll',...
        'ampVelCorrMean','ampVelCorrStd','ampVelAllCorr')
end


