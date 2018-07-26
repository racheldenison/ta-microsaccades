% getTrialsBehav.m

%% settings
experiment = 'e5'; %'e0','e3','e5','e0e3','e0e3e5','e0_nothresh',etc,'e0e3e5_nothresh'

saveFileName = sprintf('data/%s/behav_data.mat', experiment);
saveAnalysis = 1;

%% load subjects
switch experiment
    case 'e0e3'
        e0S = load(sprintf('data/%s/subjects.mat', 'e0')); 
        e3S = load(sprintf('data/%s/subjects.mat', 'e3')); 
        subjects = [e0S.subjects e3S.subjects];
    case 'e0e3e5'
        e0S = load(sprintf('data/%s/subjects.mat', 'e0'));
        e3S = load(sprintf('data/%s/subjects.mat', 'e3'));
        e5S = load(sprintf('data/%s/subjects.mat', 'e5'));
        subjects = [e0S.subjects e3S.subjects e5S.subjects];
    case 'e0e3e5_nothresh'
        e0S = load(sprintf('data/%s/subjects.mat', 'e0_nothresh'));
        e3S = load(sprintf('data/%s/subjects.mat', 'e3_nothresh'));
        e5S = load(sprintf('data/%s/subjects.mat', 'e5_nothresh'));
        subjects = [e0S.subjects e3S.subjects e5S.subjects];
    otherwise
        load(sprintf('data/%s/subjects.mat', experiment))
end

nSubjects = numel(subjects);

%% get trial accuracy for each subject
accData = [];
respCueCond = [];
for iS = 1:nSubjects
    subject = subjects{iS};
    filePaths = getSubjectDataFilePaths(upper(experiment), subject);
    
    nRuns = numel(filePaths);
    clear data
    for iRun = 1:nRuns
        data(iRun) = load(filePaths{iRun},'expt');
    end

    trials_headers = data(1).expt.trials_headers;
    cueIdx = find(strcmp(trials_headers,'cuedInterval'));
    respCueIdx = find(strcmp(trials_headers,'respInterval'));
    accIdx = find(strcmp(trials_headers,'correct'));
    fixFields = fields(data(1).expt.eye);
    nFF = numel(fixFields);

    %% concatenate run data, trials presented and fixation maintenance
    trials = [];
    fix = logical([]);
    for iRun = 1:nRuns
        trials = [trials; data(iRun).expt.trialsPresented.trials];
        ff = [];
        for iFF = 1:nFF
            ff(:,iFF) = data(iRun).expt.eye.(fixFields{iFF});
        end
        fix = [fix; all(ff,2)];
    end
    
    %% get accuracy and cue condition for all trials with no fixation break
    acc = trials(fix, accIdx);
    cueCond = trials(fix, cueIdx);
    respCue = trials(fix, respCueIdx);
    
    %% split up acc into three columns, according to cue
    cues = unique(cueCond);
    nCues = numel(cues);
    
    for iCue = 1:nCues
        cue = cues(iCue);
        
        cueByCond{iCue} = cueCond(cueCond==cue); % sanity check
        accByCond0{iCue} = acc(cueCond==cue);
        respCueByCond0{iCue} = respCue(cueCond==cue);
        nTrialsInCond(iCue) = length(accByCond0{iCue});
    end
    
    % put into matrix
    if strcmp(experiment,'e5')
        accByCond = nan(300,nCues);
        respCueByCond = nan(300,nCues);
    else
        accByCond = nan(max(nTrialsInCond),nCues);
        respCueByCond = nan(max(nTrialsInCond),nCues);
    end
    for iCue = 1:nCues
        accByCond(1:nTrialsInCond(iCue),iCue) = accByCond0{iCue};
        respCueByCond(1:nTrialsInCond(iCue),iCue) = respCueByCond0{iCue};
    end
    
    % store
    accData(:,:,iS) = accByCond;
    respCueCond(:,:,iS) = respCueByCond;
    trialsPresented{iS} = trials;
end

%% visualize accData
for i = 1:nSubjects
    imagesc(respCueCond(:,:,i))
    pause(1)
end

%% save
if saveAnalysis
    save(saveFileName, 'subjects', 'filePaths', 'trials_headers', 'cues', 'trialsPresented', 'accData', 'respCueCond')
end
