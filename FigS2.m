%% Fig S2B choose most different denovos
% ordernow from fig2.m
bigcorr = corr(StrainSumProm(Protype<3, goodids(ordernow)),'rows','pairwise');
[byrows,id] = min(bigcorr(:,:));
temp1=goodids(ordernow(1));
temp2=goodids(ordernow(19));

figure
subplot(1,2,1)
scatter(StrainSumProm(Protype<3,temp1),StrainSumProm(Protype<3,temp2),30,'k','filled','markerfacealpha',0.5)
xlabel(strains.strain(temp1))
ylabel(strains.strain(temp2))
hold on
plot(xlim,xlim,'--','Color',[.75 .75 .75])
title(num2str(corr(StrainSumProm(Protype<3,temp1),StrainSumProm(Protype<3,temp2))))

subplot(1,2,2)
temp1 = goodids(orderClust(end-3));
temp2 = goodids(orderClust(end-2));
scatter(StrainSumProm(Protype<3,temp1),StrainSumProm(Protype<3,temp2),30,'k','filled','markerfacealpha',0.5)
xlabel(strains.strain(temp1))
ylabel(strains.strain(temp2))
hold on
plot(xlim,xlim,'--','Color',[.75 .75 .75])
title(num2str(corr(StrainSumProm(Protype<3,temp1),StrainSumProm(Protype<3,temp2))))
%% Figure S2C show scatter of mean and disordered fraction
firstbatchyes = find(ismember(msn2denovo.strainid, FirstBatch.strainid(1:end-2)));
DMyes = find(ismember(msn2denovo.strainid, denovoinfo.strainid));
DMyes = DMyes(ismember(DMyes,firstbatchyes));
temp = DMyes;
figure
subplot(1,2,1)
hold on
fill([0.6 0.6 1 1],[0.9 1 1 0.9],'y')
scatter(msn2denovo.MeanPred(temp),msn2denovo.DISOF(temp),30,msn2denovo.maxcorr(temp),'filled','MarkerEdgeColor',[.5 .5 .5],'MarkerFaceAlpha',0.5)
colormap(brewermap(1000,'Purples'))
caxis([0.7 1])
colorbar
xlabel('mean predicted scores')
ylabel('fraction of disordered regions')
selected = find(msn2denovo.DISOF(temp)<0.95);
[~,idx] = sort(msn2denovo.DISOF(selected),'descend');

subplot(1,2,2)
n=1
for i =firstbatchyes(selected(idx))'
    subplot(7,1,n)
    scatter(1:numel(msn2denovo.MetaPredict{i}),msn2denovo.MetaPredict{i},10,msn2denovo.MetaPredict{i}>0.5,'filled')
         colormap([88,  89,  91; 247, 148, 29] / 255)
    ylim([0 1])
    % yline(0.5)
    tempid2 = find(ismember(denovoinfo.StrainName,msn2denovo.Sample{i}));
    tempid1 = denovoinfo.strainid(tempid2);
    title(sprintf('%s, %s, maxcorr = %.2f',msn2denovo.Sample{i}, ...
        denovoinfo.newname{tempid2},msn2denovo.maxcorr(i)))
    n=n+1;
    xlim([0 460])
end
