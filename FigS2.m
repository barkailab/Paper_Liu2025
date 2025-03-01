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