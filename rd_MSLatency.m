% rd_MSLatency.m

expDir = 'data/e5'; % check condNames
analysisName = 'bonneh_post2_dir'; %'onsets', 'bonneh_pre2'

load(sprintf('%s/analysis_struct', expDir))

if strfind(expDir, 'e5')
    condNames = {'EVENT_CUE_0', 'EVENT_CUE_1', 'EVENT_CUE_2', 'EVENT_CUE_3'};
else
    condNames = {'EVENT_CUE_0', 'EVENT_CUE_1', 'EVENT_CUE_2'};
end
nConds = numel(condNames);
eventTime = 1500;
preTime = 500; % time of cue
postTime = 2750; % end of timeseries for two-target experiments

microsaccades_data = analysis_struct.microsaccades_data;
nSubjects = numel(microsaccades_data);

for iSubject = 1:nSubjects
    onsetsAll{iSubject} = [];
    condsAll{iSubject} = [];
    for iCond = 1:nConds
        onsets = microsaccades_data{iSubject}.(condNames{iCond}).onsets;
        directions = microsaccades_data{iSubject}.(condNames{iCond}).directions;
        if iscell(onsets)
            nTrials = numel(onsets);
        else
            nTrials = size(onsets,1);
        end
        for iTrial = 1:nTrials
            if iscell(onsets)
                vals = onsets{iTrial};
            else
                vals = find(onsets(iTrial,:));
            end
            
            if isempty(vals)
                dirs = [];
            else
                dirs = directions(1:numel(vals));
                directions = directions(numel(vals)+1:end);
            end
            
            if strfind(analysisName, 'pre')
                idx = find(vals<eventTime & vals>preTime,1,'last');
                val = vals(idx); % pre
                d = dirs(idx);
            elseif strfind(analysisName, 'post')
                idx = find(vals>=eventTime & vals<postTime,1,'first');
                val = vals(idx); % post
                d = dirs(idx);
            else
                fprintf('did not find analysisName')
            end
            
            if isempty(val)
                lastOnset(iTrial,iCond,iSubject) = nan;
                lastDir(iTrial,iCond,iSubject) = nan;
            else
                lastOnset(iTrial,iCond,iSubject) = val;
                lastDir(iTrial,iCond,iSubject) = d;
            end
        end
        
        onsetsAll{iSubject} = [onsetsAll{iSubject}; onsets];
        conds = ones(size(onsets,1),1)*iCond;
        condsAll{iSubject} = [condsAll{iSubject}; conds];
    end
end

% check
if ~isempty(directions)
    error('onsets and directions may not have matched up properly!')
end

lastOnset(lastOnset==0) = nan;

lastOnsetMean = squeeze(nanmean(lastOnset,1))';

boneh_mat = lastOnsetMean;

save(sprintf('%s/%s.mat',expDir,analysisName), 'boneh_mat', 'lastOnset', 'lastDir', 'condNames', 'eventTime')
% save(sprintf('%s/%s.mat',expDir,analysisName), 'onsetsAll','condsAll','condNames')

