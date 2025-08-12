%% Fig S5 A,B
%% draw corrDBDand corrFL by each pertubation
figure
[~,order] = sort(pertubationinfo.order);% choose whatever order you want
c=0;
for i = 1:height(pertubationinfo.order)
    c=c+1;
    % corr DBD
    subplot(6,10,c-(floor((c-1)/10))*10+20*(floor((c-1)/10)))
    hold on
    cursel = find(ismember(denovoinfo.newtype,pertubationinfo.fulltype{i}));
    cursel = cursel(strains.max_corr(denovoinfo.strainid(cursel))>0.9);
    id0= find(ismember(denovoinfo.newname,[pertubationinfo.fulltype{i}(1),'F0']));

    refi = refsel(str2num(pertubationinfo.fulltype{i}(1)));

    % get x axis variable
    try
        xn = denovoinfo.(pertubationinfo.var{i})(cursel);
        usdedvar = pertubationinfo.var{i};
    catch
        xn = denovoinfo.nST(cursel);
        usdedvar = 'nST';
    end
    xn = denovoinfo.(usdedvar)(refi)-xn;

    DBDtemp = drawlist(15,2);
    scatter(xn,corr(StrainSumProm(Protype<3, ...
        denovoinfo.strainid(cursel)),StrainSumProm(Protype<3,DBDtemp)), ...
        5+xn/2,denovoinfo.mutationid(cursel),'filled')
    scatter(0,corr(StrainSumProm(Protype<3, ...
        denovoinfo.strainid(refi)),StrainSumProm(Protype<3,DBDtemp)), ...
        5,18,'filled')
    plot([0;xn],corr(StrainSumProm(Protype<3, ...
        denovoinfo.strainid([refi;cursel])),StrainSumProm(Protype<3,DBDtemp)), ...
        'color',allcolors(denovoinfo.mutationid(cursel(1)),:))
    colormap(gca,allcolors)
    caxis([1 22])
    title(pertubationinfo.fulltype{i})
    yline(0.9,'--','Color',[.5 .5 .5])
    set(gca,'XDir','normal')
    ylim([ 0.5 0.75])
    % corr Reference
    subplot(6,10,c-10*(floor((c-1)/10))+20*(floor((c-1)/10))+10)
    hold on

    scatter(xn,corr(StrainSumProm(Protype<3, ...
        denovoinfo.strainid(cursel)),StrainSumProm(Protype<3,denovoinfo.strainid(refi))), ...
        5+xn/2,denovoinfo.mutationid(cursel),'filled')

    scatter(0,corr(StrainSumProm(Protype<3, ...
        denovoinfo.strainid(refi)),StrainSumProm(Protype<3,denovoinfo.strainid(refi))), ...
        5,18,'filled')
    plot([0;xn],corr(StrainSumProm(Protype<3, ...
        denovoinfo.strainid([refi;cursel])),StrainSumProm(Protype<3,denovoinfo.strainid(refi))), ...
        'color',allcolors(denovoinfo.mutationid(cursel(1)),:))
    colormap(gca,allcolors)
    caxis([1 22])
    yline(0.9,'--','Color',[.5 .5 .5])
    xlabel(['n',pertubationinfo.type{i}])
    set(gca,'XDir','normal')
    ylim([0.8 1])
end
exportgraphics(gcf,'Fig6AB.pdf','ContentType','vector')
%% FigS5C abudnace and corrDBD corrFL
figure
hold on
tempsel = find(denovoinfo.abd~=0);
scatter(denovoinfo.xRef(tempsel),denovoinfo.abd(tempsel),40,denovoinfo.xDBD(tempsel),'filled')
[x,y] = polyfit(denovoinfo.xRef(tempsel),denovoinfo.abd(tempsel),1);

x_fit = linspace(0,1,100);
y_fit = polyval(x, x_fit);
plot(x_fit, y_fit, '-', 'color',[.5 .5 .5],'LineWidth', 2);

colormap(brewermap(1000,'Blues'))
caxis([0 1])
colorbar
xlabel('Corr. Reference')
ylabel('Abandance')
title('c = 0')
%% Fig S5E compare de novo motifs of Denovo-DBD with Msn2
tempsel = msn2denovo.strainid(msn2denovo.maxcorr>=0.9);
DDselinDK = find(ismember([DK_normP.strainid], tempsel));

mybinedges = 0:0.02:1;
figure 
histogram(UnbMCI(DDselinDK),'Normalization','probability','BinEdges',mybinedges,'FaceColor','#349946','FaceAlpha',0.6,'EdgeColor','none')
hold on
histogram(cat(2,withother{DDselinDK}),'Normalization','probability','BinEdges',mybinedges,'FaceColor','#a8cedf','FaceAlpha',1,'EdgeColor','none')
legend({'to in-vitro Msn2 motif','to other in-vitro motifs'})
ylabel("Probability")
title('Denovo-DBD to other in vitro')
xlabel('Euclidean motif distance')
%% Fig S5F fold change of absolute promoter measurement and fraction of motifs-containing targets
figure
imagesc(PerbMat,[0.3 1])
set(gca,'ytick',1:height(pertubationinfo),'yticklabel',pertubationinfo.fulltype)
colormap([.5 .5 .5;brewermap(1000,'Greens')])
hold on
xline(0.5:1:4.5,'white')
yline(0.5:1:24.5,'white')
title('Fraction of motif appearance in targets')
ylim([0.5 21.5])
xticks(1:5)
figure
imagesc(PerbMatMotifFra,[0.3 1])
set(gca,'ytick',1:height(pertubationinfo),'yticklabel',pertubationinfo.fulltype)
colormap([.5 .5 .5;brewermap(1000,'Greens')])
hold on
xline(0.5:1:4.5,'white')
yline(0.5:1:24.5,'white')
title('Fraction of motif appearance in targets')
ylim([0.5 21.5])
xticks(1:5)
%% Fig S5G scatter of fraction of motifs-containing targets and correlation with reference
figure
scatter(denovoinfo.MotifFra(sel), denovoinfo.xRef(sel),50,denovoinfo.Pratio(sel),'filled','MarkerEdgeColor',[.5 .5 .5],'MarkerFaceAlpha',0.5)
colormap(brewermap(100,'-RdBu'))
clim([-0.5 0.5])
xlabel('Fraction of motif-containing promoters')
ylabel('Corr. Ref')
c=colorbar;
c.Label.String = 'log(motif binding./Ref motif binding)';
%% Fig S5H number of targets bound by native TF, de novo-DBD
cur_sp = StrainSumProm;
thr=3;
log_sum = log2(cur_sp+700); % supress noise
zscore= (log_sum-mean(log_sum,1,'omitnan'))./std(log_sum,[],1,'omitnan');
selMD = (msn2denovo.strainid(msn2denovo.maxcorr>=0.9));
selDD = acol(drawlist(DDid([1:6,8,9]),31:36));
anyMDsig = any(zscore(:,selMD)>thr&Protype<3,2);
anyDDsig = any(zscore(:,[selMD;selDD])>thr&Protype<3,2);
zMsn2 = zscore(:,281);
Msn2sig = (zMsn2>thr&Protype<3);

zallTF =  (TFsumProm-mean(TFsumProm,1,'omitnan'))./std(TFsumProm,[],1,'omitnan');
z_TF_sig = zallTF>thr&Protype<3;
selTFsig = (sum(z_TF_sig,2)>2);

DD0TFsig = selTFsig.*anyDDsig;
MD0TFsig = selTFsig.*anyMDsig;
MD0DDsig = anyDDsig.*anyMDsig;
%% plot
figure
data = [sum(DD0TFsig), sum(selTFsig)-sum(DD0TFsig); sum(DD0TFsig), sum(anyDDsig)-sum(DD0TFsig);
    sum(MD0TFsig),sum(Msn2sig)-sum(MD0TFsig);sum(MD0TFsig),sum(anyMDsig)-sum(MD0TFsig)];
data2 = [0, sum(selTFsig); sum(DD0TFsig), sum(anyDDsig)-sum(DD0TFsig);
    sum(MD0TFsig),sum(anyMDsig)-sum(MD0TFsig)];
% data2 = [data2,[0;0;]]
bar(data2,'stacked','EdgeColor','none')
legend({' ','Common targets with native TF'})
xticklabels({'Native TF targets','All de novo-DBD',"All de novo-DBD^Msn2"})
%% Fig S5I predicted disorder tendency
firstbatchyes = find(ismember(msn2denovo.strainid, FirstBatch.strainid(1:end-2)));
DMyes = find(ismember(msn2denovo.strainid, denovoinfo.strainid));
DMyes = DMyes(~ismember(DMyes,firstbatchyes));
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
