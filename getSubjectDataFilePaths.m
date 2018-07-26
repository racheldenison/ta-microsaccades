function paths = getSubjectDataFilePaths(study, subjectInit)

% from jake's repo: https://github.com/jparkswims/tempatten
% tempatten2.m
% - ta_params.m
% - ta_preprocess.m
% -- trialsx.m

%% setup
% study = 'E0';
% subjectInit = [];

if nargin==0
    error('must input study as first arg: E0, E3, E5')
end
if nargin<2
    subjectInit = [];
end

pathOption = 'fullpath'; % 'filename','fullpath'

%% params
switch study
    case 'E0'
        %     subjects = {'ma' 'ad' 'bl' 'ec' 'ty' 'zw' 'hl' 'rd' 'jp'}; % Pupil, no VP
        %     runs = [2 2 2 2 2 2 4 4 4];
        subjects = {'ma' 'ad' 'bl' 'ec' 'ty' 'vp' 'zw' 'hl' 'rd' 'jp'}; % MS, includes VP
        runs = [2 2 2 2 2 2 2 4 4 4];
        locs = [0 1000 1250 1750];
%         trials = 640;
    case 'E3'
        subjects = {'bl' 'ca' 'ec' 'en' 'ew' 'id' 'jl' 'jx' 'ld' 'ml' 'rd' 'sj'}; % Note, id not used for MS
        runs = [4 4 4 4 4 4 4 4 4 4 4 4];
        locs = [0 1000 1250 1750];
%         trials = 640;
    case'E5'
        subjects = {'ds' 'gb' 'gb2' 'ht' 'ik' 'jg' 'jp' 'rd' 'xw' 'yz'}; % Note, gb2 not used for MS
        runs = [1 1 1 1 1 1 1 1 1 1];
        locs = [0 1000 1250 1500 2000];
%         trials = 960;
    otherwise
        error('study not recognized')
end

%% paths
if strcmp(study,'E0')
%     TAeyepaths = {'/Volumes/purplab/EXPERIMENTS/1_Current Experiments/Rachel/Temporal_Attention/eyedata/E0_cb/'...
%         '/Volumes/purplab/EXPERIMENTS/1_Current Experiments/Rachel/Temporal_Attention/eyedata/E2_SOA_cbD6/'...
%         '/Volumes/purplab/EXPERIMENTS/1_Current Experiments/Rachel/Temporal_Attention/eyedata/pilot/'};
    TAdatapaths = {'/Volumes/purplab/EXPERIMENTS/1_Current Experiments/Rachel/Temporal_Attention/data/E0_cb/'...
        '/Volumes/purplab/EXPERIMENTS/1_Current Experiments/Rachel/Temporal_Attention/data/E2_SOA_cbD6/'...
        '/Volumes/purplab/EXPERIMENTS/1_Current Experiments/Rachel/Temporal_Attention/data/pilot/'};
%     pathsub = {[1:6] [7 8] [9]}; % Pupil
    pathsub = {[1:7] [8 9] [10]}; % MS
%     edffind = '%s/%s/*%d_run0%d*.edf';
    matfind = '%s/%s/*%d_run0%d*Temp*.mat';
elseif strcmp(study, 'E3')
%     TAeyepaths = {'/Volumes/purplab/EXPERIMENTS/1_Current Experiments/Rachel/Temporal_Attention/eyedata/E3_adjust/'};
    TAdatapaths = {'/Volumes/purplab/EXPERIMENTS/1_Current Experiments/Rachel/Temporal_Attention/data/E3_adjust/'};
    pathsub = {[1:12]};
%     edffind = '%s/%s/*%d_run0%d*.edf';
    matfind = '%s/%s/*%d_run0%d*Temp*.mat';
elseif strcmp(study, 'E5')
%     TAeyepaths = {'/Volumes/purplab/EXPERIMENTS/1_Current Experiments/Rachel/Temporal_Attention/eyedata/E5_T3_cbD15/'};
    TAdatapaths = {'/Volumes/purplab/EXPERIMENTS/1_Current Experiments/Rachel/Temporal_Attention/data/E5_T3_cbD15/'};
    pathsub = {[1:10]};
%     edffind = 'lol';
    matfind = '%s/%s/*run01WW*Temp*.mat';
end

%% get data file paths for each subject and run
matfiles = {};
for p = 1:length(TAdatapaths)
    for s = pathsub{p}
%         TAeyepath = TAeyepaths{p};
        TAdatapath = TAdatapaths{p};
        % pa = trialsx(pa,TAeyepath,TAdatapath,edffind,matfind,pa.subjects{s},study,type,s);
        subject = subjects{s};
        
        for r = 1:runs(s)
            if strcmp(study,'E0') || strcmp(study,'E3') || strcmp(study,'E0E3')
                b = dir(sprintf(matfind,TAdatapath,subject,locs(end-1),r));
            elseif strcmp(study,'E5')
                b = dir(sprintf(matfind,TAdatapath, subject));
            end
            
            if length(b) == 1
                matfile = b.name;
            else
                matfile = b(length(b)).name;
            end
            
            matfiles{s,r} = matfile;
            matfilesfull{s,r} = sprintf('%s%s/%s', TAdatapath, subject, matfile);
        end
    end
end

%% return path for requested subject
switch pathOption
    case 'filename'
        pp = matfiles;
    case 'fullpath'
        pp = matfilesfull;
    otherwise
        error('pathOption not recognized')
end

if isempty(subjectInit)
    paths = pp;
else
    s = strcmp(subjectInit, subjects);
    paths = pp(s,1:runs(s));
    
    if ~any(s)
        warning('subject initials did not match any subject in this study')
    end
end

