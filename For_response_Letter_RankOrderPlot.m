%% rank order plot
 
[sortedSignal,idx]=sort(StrainSumProm(Protype<3,1),1,'descend');
log_sum = log2(StrainSumProm+700); 
zscore= (log_sum-mean(log_sum,1,'omitnan'))./std(log_sum,[],1,'omitnan');
zscore = zscore(Protype<3);
% sortedSignal = StrainSumProm(Protype<3,:);
sortedMotif = motifYesMat(idx,1);
cuttoff = max(find(zscore(idx)>3));
figure
hold on
fill([0 0 cuttoff cuttoff], [0 18 18 0],'y','edgecolor','none')
scatter(1:5358,log2(sortedSignal),10,motifYesMat(idx),'filled')
colormap(brewermap(2,'-Spectral'))
legend('targets','Motif+', 'Motif-')
xlabel('Rank'); ylabel('ChEC signal'); title('Ranked ChEC signal by motif')
title('MSN2_aG89_aP114_aN43_aQ53_aF63.mat')
