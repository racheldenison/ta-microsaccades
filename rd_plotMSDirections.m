% rd_plotMSDirections.m

%% load data
experiment = 'e3'; %'e0','e3'

switch experiment
    case 'e0'
        load data/e0/analysis_struct4.mat
        subjectNames = {'ad','bl','ec','hl','jp','ma','rd','ty','vp','zw'};
    case 'e3'
        load data/e3/analysis_struct-11_subjects.mat
        subjectNames = {'o1','o2','o3','o4','o5','o6','o7','o8','o9','o10','o11'};
    otherwise
        error('experiment not found')
end

% {bin, condition, subject}

%% setup
timeBins = {-999:-800, -799:-600, -599:-400, -399:-200, -199:0, 1:200, 201:400, 401:600, 601:800};
condNames = {'attend both','attend t1','attend t2'};

nDirBins = 8;
dirBinSize = 2*pi/nDirBins;
dirBinEdges = (-dirBinSize:dirBinSize:2*pi-dirBinSize) + dirBinSize/2;

directions = analysis_struct.results_grand_total.directions;

[nTimeBins, nConds, nSubjects] = size(directions);
binDuration = timeBins{1}(end)-timeBins{1}(1)+1;

nTrials = [.2 .4 .4]*640;

%% input checks
if nTimeBins~=numel(timeBins)
    error('check number of timeBins')
end
if nConds~=numel(condNames)
    error('check number of condNames')
end
if nSubjects~=numel(subjectNames)
    error('check number of subjects')
end
nTrialsCheck = [numel(analysis_struct.microsaccades_data{1}.EVENT_CUE_0.number_of_saccades), ...
    numel(analysis_struct.microsaccades_data{1}.EVENT_CUE_1.number_of_saccades), ...
    numel(analysis_struct.microsaccades_data{1}.EVENT_CUE_2.number_of_saccades)];
if any(nTrialsCheck~=nTrials)
    error('check number of trials per condition')
end

%% count directions
counts = zeros(length(dirBinEdges)-1, nTimeBins, nConds, nSubjects);

for iSubject = 1:nSubjects
    for iCond = 1:nConds
        for iTimeBin = 1:nTimeBins
            d = directions{iTimeBin,iCond,iSubject};
            h = polarhistogram(d,dirBinEdges);
            counts(:,iTimeBin,iCond,iSubject) = h.Values;
        end
    end
end

%% convert to rates
rates = zeros(size(counts));
for iCond = 1:nConds
    rates(:,:,iCond,:) = counts(:,:,iCond,:)/(binDuration/1000*nTrials(iCond));
end

%% average counts or rates per direction
groupData = rates;
groupMean = mean(groupData,4);
groupCondMean = mean(groupMean,3);

%% plot mean across conditions
rlims = [0 .7]; % [0 .3]
figure('Position',[150 650 2000 550])
for iTimeBin = 1:nTimeBins
    subplot(1,nTimeBins,iTimeBin)
    h = polarhistogram('BinEdges',dirBinEdges,'BinCounts',groupCondMean(:,iTimeBin));
    h.Parent.RLim = rlims;
    title(sprintf('%d-%d ms',timeBins{iTimeBin}(1),timeBins{iTimeBin}(end)))
end

%% plot mean per condition
rlims = [0 .8]; %[0 .4];
figure('Position',[150 650 2000 550])
for iTimeBin = 1:nTimeBins
    subplot(1,nTimeBins+1,iTimeBin)
    for iCond = 1:nConds
        h = polarhistogram('BinEdges',dirBinEdges,'BinCounts',groupMean(:,iTimeBin,iCond));
        h.DisplayStyle = 'stairs';
        h.LineWidth = 2;
        h.Parent.RLim = rlims;
        hold on
    end
    title(sprintf('%d-%d ms',timeBins{iTimeBin}(1),timeBins{iTimeBin}(end)))
end
% plot legend
subplot(1,nTimeBins+1,iTimeBin+1)
for iCond = 1:nConds
    h = polarhistogram('BinEdges',dirBinEdges,'BinCounts',ones(1,length(dirBinEdges)-1));
    h.DisplayStyle = 'stairs';
    h.LineWidth = 2;
    h.Parent.RLim = rlims;
    hold on
end
title('[Hz]')
legend(condNames)

%% plot mean per subject
rlims = [0 1];
figure('Position',[400 50 800 1300])
for iSubject = 1:nSubjects
    for iTimeBin = 1:nTimeBins
        subplot(nSubjects,nTimeBins+1,iTimeBin+(nTimeBins+1)*(iSubject-1))
        for iCond = 1:nConds
            h = polarhistogram('BinEdges',dirBinEdges,'BinCounts',groupData(:,iTimeBin,iCond,iSubject));
            h.DisplayStyle = 'stairs';
            h.LineWidth = 2;
            h.Parent.RLim = rlims;
            hold on
        end
        if iSubject==1
            title(sprintf('%d-%d ms',timeBins{iTimeBin}(1),timeBins{iTimeBin}(end)))
        end
        if iSubject~=1 || iTimeBin~=1
            h.Parent.ThetaTickLabel = {};
            h.Parent.RTickLabel = {};
        end
    end
    subplot(nSubjects,nTimeBins+1,iTimeBin+1+(nTimeBins+1)*(iSubject-1))
    title(subjectNames{iSubject})
    axis off
end


