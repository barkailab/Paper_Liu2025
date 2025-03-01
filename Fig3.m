DD_sel = [21,22,23,19,2,24]; % from FirstBatch
DDid = find(drawlist(:,31)>0);
%% Fig3B draw corr. DBD and WT
shiftScore = (corrWT(DDid(1:9),31:36)-corrDBD(DDid(1:9),31:36))./(2-corrWT(DDid(1:9),2));

figure
hold on
% hold on
for i = 1:9
    hold on
    scatter(i*ones(6,1),1:6,20+100*adjmaxcorr(DDid(i),31:36),shiftScore((i),:),'filled','MarkerEdgeColor','#635A71')
    colormap(brewermap(1000,'PuOr'))

    caxis([-0.3 0.3])
    xlim([0 11])
    ylim([0 7])
    set(gca,'YDir','reverse')
    xticks([])
    yticks([])
end
%% write to table S5 about values
temptable = table();
temptable.denovo = strrep(strrep(denovolist_new,'a',' '),'_','');
for i = 1:9
    temptable.([TFinfo.TF{DDid(i)},'xDBD-only']) = corrDBD(DDid(i),31:36)';
    temptable.([TFinfo.TF{DDid(i)},'xFL']) = corrWT(DDid(i),31:36)';
    temptable.([TFinfo.TF{DDid(i)},' Shift Score']) = shiftScore((i),:)';

end

writetable(temptable,'/home/labs/barkailab/jingliu/Documents/Finalfiles/Sfiles/SuppValues.xlsx','Sheet','Figure3B')
%% Fig3B first row DBDxFL
ax = axes;
scatter(ax,1:9,zeros(9,1)*7,10+50,corrDBD(DDid(1:9),1),'filled','MarkerEdgeColor','#635A71')
colormap(ax,brewermap(1000,'Purples'))
ylim([0 7])
xlim(ax,[0 11])

    set(gca,'YDir','reverse')

xticks(1:9)
xticklabels(TFinfo.TF(DDid(1:9)))
%% Fig 3 C corr bindscore & reproducibility
figure
hold on
for i = 1:9
    scatter(maxcorr(denovoTFid(i),31:36),(DDbindScore(i,3:8)./DDbindScore(i,1)),30,corrWT(denovoTFid(i),31:36),'filled','MarkerEdgeColor',[0.2471         0    0.4902])
    colormap(brewermap(1000,'Purples'))
    caxis([0 1])
end
xlabel('r: repeats')
ylabel('In vitro motifs binding score./WT')
title('Reproducibility & binding in vitro motifs')
xline(0.9,'--')
colorbar

%% select reproducible ones
fulllist = acol(drawlist(DDid(1:9),[2,31:36]));
goodsel = acol(maxcorr(DDid(1:9),[2,31:36]))>0.9;

allTFlabels = acol(repmat(TF(DDid(1:9)),6,1)');
%% Fig 3D
figure
allxDBD = acol(corrDBD(DDid(1:9),[31:36]));
allxFL = acol(corrWT(DDid(1:9),[31:36]));
goodsel = acol(maxcorr(DDid(1:9),[31:36]))>0.9;

scatter(allxDBD(goodsel),allxFL(goodsel),60,'blue','filled','MarkerFaceAlpha',0.5)
hold on
plot(xlim,xlim,'--','Color',[.75 .75 .75])
xlabel('Corr. DBD')
ylabel('Corr. FL')
title('Fig3E')
%% Fig 3E big heatmap corr by each DBD
figure
imagesc(corr(StrainSumProm(Protype<3,fulllist(goodsel)),'rows','pairwise'),[0,1])
set(gca,'ytick',movmean([1,goodseln],2,'Endpoints','discard'),'YTickLabel',[{'DBD'},denovolist_new'],'xtick',[]);
colormap(gca,brighten(flipud(brewermap(1000,'-Blues')),0.3))
yline(goodseln+0.5,'k')
xline(goodseln+0.5,'k')
title('corr DBD-Denovo')
colorbar
%% Fig3G WT corr nd reproducibility
figure
scatter(ones(6,1),1:6,10+50*adjmaxcorr(27,31:36),corrWT(27,31:36),'filled')
colormap(brewermap(1000,'Purples'))
caxis([0 1])
xlim([0.5 2.5])
ylim([0 7])
set(gca,'YDir','reverse')
xticks([])
yticks([])