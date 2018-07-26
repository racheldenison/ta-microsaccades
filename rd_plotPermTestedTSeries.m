% rd_plotPermTestedTSeries.m

%% setup
filename = '/Users/rachel/Google Drive/NYU/Projects/Temporal_Attention/Code/Microsaccade/data/combined/TimelinesPermutationTests_E0E3E5_msrate_-500to0ms_N30_nothresh_perm0.txt';

%% read data and format headers
data = dlmread(filename,' ',1,1);

fileID = fopen(filename);
C = textscan(fileID,'%s',size(data,2));
fclose(fileID);

headers = C{1};
for i = 1:numel(headers)
    headers{i} = headers{i}(2:end-1);
end
betaHeaders = headers(3:10);

%% clusters
clusters.beta_cueneutral = dlmread('data/combined/TimelinesPermutationTests_E0E3E5_msrate_-500to0ms_N30_nothresh_beta_cueneutral_noheader.txt');

%% time
t = data(:,strcmp(headers,'time'));

%% significant beta values
for iBeta = 1:numel(betaHeaders)
    betaName = betaHeaders{iBeta};
    betaData(:,iBeta) = data(:,strcmp(headers,betaName));
    betaDataSig(:,iBeta) = logical(data(:,strcmp(headers,sprintf('%s.significant.permtested.neg', betaName))) ...
        + data(:,strcmp(headers,sprintf('%s.significant.permtested.pos', betaName))));
end

%% plot significant beta and interaction timeseries
betaName = 'beta_cueneutral'; % 'beta_cueneutral', 'beta_cuet2'

betaVals = betaData(:,strcmp(betaHeaders,betaName));
betaSig = betaDataSig(:,strcmp(betaHeaders,betaName));
bs = nan(size(betaVals));
bs(betaSig) = betaVals(betaSig);
betaNameE3 = sprintf('%s_e3', betaName);
betaNameE5 = sprintf('%s_e5', betaName);
e3Sig = betaDataSig(:,strcmp(betaHeaders,betaNameE3));
e5Sig = betaDataSig(:,strcmp(betaHeaders,betaNameE5));

%% plot
clusterColors = {[.8 .8 .8]};
plotClusters = 1;
ylims = [-0.2 1];

figure
hold on
if plotClusters
    vals = clusters.(betaName);
    for iC = 1:size(vals,1)
        clusterTimes = vals(iC,2:3);
        rectangle('Position',[clusterTimes(1), ylims(1), diff(clusterTimes), diff(ylims)],...
            'FaceColor', clusterColors{1}, 'EdgeColor', 'none')
    end
end

plot(t, betaVals)
plot(t, bs,'LineWidth',4)
plot(t, e3Sig)
plot(t, e5Sig)
xlabel('time (ms)')
ylabel('beta')
legend({und2space(betaName), 'significant', 'e3 interaction', 'e5 interaction'})
