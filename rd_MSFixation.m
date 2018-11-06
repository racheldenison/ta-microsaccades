% rd_MSFixation.m

%% settings
experiment = 'e0e3e5'; %'e0','e3','e5','e0e3e5'

saveFileName = sprintf('data/%s/FixationAnalysis.mat', experiment);
saveAnalysis = false;

%% load data
switch experiment
    case 'e0'
        load data/e0/analysis_struct_xy_1.mat
        load(sprintf('data/%s/subjects.mat', experiment))
        screen = load(sprintf('data/%s/screen.mat', experiment)); % from rd_analyzeBehavPupilMSSubjects
        exp = repmat(cellstr(''),1,numel(subjects));
        t = -1000:1749; % where 0 is time of T1 (Omer segmented these trials starting at time of precue)
        twin = [-1000 750-120]; % [precue response cue-120 ms eye slack]
        sz = [640 2750]; 
    case 'e3'
        load data/e3/analysis_struct_xy_2.mat
        load(sprintf('data/%s/subjects.mat', experiment)) 
        screen = load(sprintf('data/%s/screen.mat', experiment));
        exp = repmat(cellstr('E3'),1,numel(subjects));
        t = -1000:1749;
        twin = [-1000 750-120];
        sz = [640 2750]; 
    case 'e5'
        load('data/e5/analysis_struct_xy_3.mat');
        load(sprintf('data/%s/subjects.mat', experiment))
        screen = load(sprintf('data/%s/screen.mat', experiment));
        exp = repmat(cellstr('E5'),1,numel(subjects));
        t = -1000:3499;
        twin = [-1000 1000]; % fixation was required until the go cue, so let's end at response cue exactly
        sz = [960 4500]; 
    case 'e0e3e5'
        e0 = load('data/e0/analysis_struct_xy_1.mat');
        e3 = load('data/e3/analysis_struct_xy_2.mat');
        e5 = load('data/e5/analysis_struct_xy_3.mat');
        analysis_struct.microsaccades_data = ...
            [e0.analysis_struct.microsaccades_data, e3.analysis_struct.microsaccades_data, ...
            e5.analysis_struct.microsaccades_data];
        e0S = load(sprintf('data/%s/subjects.mat', 'e0'));
        e3S = load(sprintf('data/%s/subjects.mat', 'e3'));
        e5S = load(sprintf('data/%s/subjects.mat', 'e5'));
        subjects = [e0S.subjects e3S.subjects e5S.subjects];
        e0Scr = load(sprintf('data/%s/screen.mat', 'e0'));
        e3Scr = load(sprintf('data/%s/screen.mat', 'e3'));
        e5Scr = load(sprintf('data/%s/screen.mat', 'e5'));
        fnames = fieldnames(e3Scr.screenInfo);
        for iF = 1:numel(fnames)
            f = fnames{iF};
            screen.screenInfo.(f) = [e3Scr.screenInfo.(f); e5Scr.screenInfo.(f)];
%             screen.screenInfo.(f) = [e0Scr.screenInfo.(f); e3Scr.screenInfo.(f); e5Scr.screenInfo.(f)];
        end
        screen.subjectInits = [e0Scr.subjectInits e3Scr.subjectInits e5Scr.subjectInits];
        exp = [repmat(cellstr(''),1,numel(e0S.subjects)) ...
            repmat(cellstr('E3'),1,numel(e3S.subjects)) ...
            repmat(cellstr('E5'),1,numel(e5S.subjects))];
        % only take the window that all experiments have
        t = -1000:1749;
        twin = [-1000 750-120];
        sz = [960 2750]; 
    otherwise
        load(sprintf('data/%s/analysis_struct.mat', experiment))
        load(sprintf('data/%s/subjects.mat', experiment))
end

%% setup
measure = 'raw_eye_data';

dataStruct = analysis_struct.microsaccades_data;
nSubjects = numel(dataStruct);

%% get screen info for each subject
for iSubject = 1:nSubjects
    subject = subjects{iSubject};
    idx = find(strcmp(screen.subjectInits, sprintf('%s%s',subject,exp{iSubject})));
    screenSize(iSubject,:) = screen.screenInfo.size(idx,:);
    screenRes(iSubject,:) = screen.screenInfo.res(idx,:);
    viewDist(iSubject,1) = screen.screenInfo.viewDist(idx,:);
    pixelsPerDegree(iSubject,1) = ang2pix(1, screenSize(iSubject,1), screenRes(iSubject,1), viewDist(iSubject), 'central');
end

cx = screenRes(:,1)/2; 
cy = screenRes(:,2)/2;

%% organize data
data = [];
for iSubject = 1:nSubjects
    conds = fields(dataStruct{iSubject}); % different conds for 2-target and 3-target subjects
    nConds = numel(conds);
    
    for iCond = 1:nConds
        rawData = dataStruct{iSubject}.(conds{iCond}).(measure);
        nTrials = numel(rawData);
        if rawData(1).x.right_eye(1)==1
            eyeName = 'left_eye';
        else
            eyeName = 'right_eye';
        end
        
        for iTrial = 1:nTrials
            data{iCond,iSubject}.(measure).x(iTrial,:) = rawData(iTrial).x.(eyeName);
            data{iCond,iSubject}.(measure).y(iTrial,:) = rawData(iTrial).y.(eyeName);
        end
    end
    eyeNames{iSubject} = eyeName;
end

%% aggregate data by subject
for iSubject = 1:nSubjects
    conds = fields(dataStruct{iSubject}); 
    nConds = numel(conds);
    
    groupData.(measure){iSubject}.x = [];
    groupData.(measure){iSubject}.y = [];
    for iCond = 1:nConds
        groupData.(measure){iSubject}.x = [groupData.(measure){iSubject}.x; data{iCond,iSubject}.(measure).x];
        groupData.(measure){iSubject}.y = [groupData.(measure){iSubject}.y; data{iCond,iSubject}.(measure).y];
    end
end

%% aggregate into trial x time x subject matrices
% E5 has a few subjects missing a couple of trials, so let's initialize
% with nan. Also, let's select position data from 1:sz(2) so that we can
% combine across experiments with same size data.
groupData.x = nan([sz nSubjects]);
groupData.y = nan([sz nSubjects]);
for iSubject = 1:nSubjects
    nTrials = size(groupData.(measure){iSubject}.x,1);
    groupData.x(1:nTrials,:,iSubject) = groupData.(measure){iSubject}.x(:,1:sz(2));
    groupData.y(1:nTrials,:,iSubject) = groupData.(measure){iSubject}.y(:,1:sz(2));
end

%% calculate distance from center of screen
cxshift = shiftdim(cx,-2);
cxmat = repmat(cxshift,sz(1),sz(2));
cyshift = shiftdim(cy,-2);
cymat = repmat(cyshift,sz(1),sz(2));
ppdshift = shiftdim(pixelsPerDegree,-2);
ppdmat = repmat(ppdshift,sz(1),sz(2));

% with reference to screen center
groupData.xpix = groupData.x-cxmat;
groupData.ypix = groupData.y-cymat;
groupData.xdeg = groupData.xpix./ppdmat;
groupData.ydeg = groupData.ypix./ppdmat;

groupData.rpix = sqrt((groupData.xpix).^2 + (groupData.ypix).^2);
groupData.rdeg = groupData.rpix./ppdmat;

%% calculate std per trial in x and y
tidx(1) = find(t==twin(1));
tidx(2) = find(t==twin(2));
xdata = groupData.xdeg(:,tidx(1):tidx(2),:);
ydata = groupData.ydeg(:,tidx(1):tidx(2),:);

% average within-trial standard deviation
groupData.xdegTrialStd = nanmean(squeeze(std(xdata,0,2)));
groupData.ydegTrialStd = nanmean(squeeze(std(ydata,0,2)));

for iSubject = 1:nSubjects
    % std of positions aggregated across all trials, all time points
    xx = xdata(:,:,iSubject);
    yy = ydata(:,:,iSubject);
    groupData.xdegExptStd(iSubject) = nanstd(xx(:));
    groupData.ydegExptStd(iSubject) = nanstd(yy(:));
    
    % within-trial mean
    groupData.xdegMean(:,iSubject) = squeeze(mean(groupData.xdeg(:,tidx(1):tidx(2),iSubject),2));
    groupData.ydegMean(:,iSubject) = squeeze(mean(groupData.ydeg(:,tidx(1):tidx(2),iSubject),2));
end

% std of the within-trial means
groupData.xdegStdOfMean = nanstd(groupData.xdegMean);
groupData.ydegStdOfMean = nanstd(groupData.ydegMean);

%% plot trials
figure
for iSubject = 1:nSubjects
    cla
    imagesc(groupData.xdeg(:,:,iSubject),[-2.5 2.5])
    title(num2str(iSubject))
    pause(2)
end

%% plot x,y
thresh = 1.5;

figure
hold on
for iSubject = 1:nSubjects
%     figure
    plot(groupData.xdegMean(:,iSubject), groupData.ydegMean(:,iSubject), '.')
    xlim([-2.5 2.5])
    ylim([-2.5 2.5])    
end
x = linspace(-thresh,thresh,100);
y = sqrt(thresh^2 - x.^2);
plot(x, y, '--k')
plot(x, -y, '--k')
axis square
xlabel('horizontal position (deg)')
ylabel('vertical position (deg)')

%% plot distance
figure
for iSubject = 1:nSubjects
    cla % for E5, too much for graphics to handle if all figs are plotted
    hold on
    plot(t, groupData.rdeg(:,:,iSubject)')
    plot(twin,[thresh thresh], 'k','LineWidth',2)
    %     xlim([0 1])
    ylim([0 3])
    title(num2str(iSubject))
    pause(2) 
end
xlabel('time (ms)')
ylabel('distance from screen center (deg)')

%% plot proportion of trials exceeding some threshold distance
vals = groupData.rdeg>thresh;

figure
hold on
plot(twin,[.5 .5], 'k','LineWidth',2)
plot(t, squeeze(mean(vals)))
xlabel('time (ms)')
ylabel('proportion exceeding threshold')



