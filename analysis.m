for iSub=1:10
vec_event_cue_0=analysis_struct.results_per_subject{iSub}.microsaccadic_rate.EVENT_CUE_0;
vec_event_cue_1=analysis_struct.results_per_subject{iSub}.microsaccadic_rate.EVENT_CUE_1;
vec_event_cue_2=analysis_struct.results_per_subject{iSub}.microsaccadic_rate.EVENT_CUE_2;
mat_0(iSub,:)=vec_event_cue_0;
mat_1(iSub,:)=vec_event_cue_1;
mat_2(iSub,:)=vec_event_cue_2;
end

% tidx = 600:1000;
tidx = 1250:1500;
for iSub=1:10
vals(iSub,1)=mean(mat_0(iSub,tidx)); 
vals(iSub,2)=mean(mat_1(iSub,tidx)); 
vals(iSub,3)=mean(mat_2(iSub,tidx)); 
end


timeAxis=[-1000:699];
figure
plot(timeAxis, mean(mat_1),timeAxis, mean(mat_2),timeAxis, mean(mat_0)); 
% zerolines
legend('attend T1','attend T2','neutral')
hold on
yzero = plot([ 0 0], get(gca,'ylim'),'k:'); % yline
yzero = plot([ 250 250], get(gca,'ylim'),'k:'); % yline
xlabel('time (ms)')
ylabel('microsaccade rate')

valsMean = mean(vals,1);
valsSte = std(vals,0,1)./sqrt(size(vals,1));

valsMean = mean(vals,1);
valsSte = std(vals,0,1)./sqrt(size(vals,1));

figure
errorbar(valsMean([2 1 3]), valsSte([2 1 3]), '.')
set(gca,'XTick',1:3)
set(gca,'XTickLabel',{'att T1','neutral', 'att T2'})

figure
plot(vals(:,[2 1 3])')
xlim([.5 3.5])
set(gca,'XTick',1:3)
set(gca,'XTickLabel',{'att T1','neutral', 'att T2'})

