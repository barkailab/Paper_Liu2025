%% A 
% mark motifs
BStable = BStables.MSN2;
sel = (BStable.type==1);
MotifsApp=BStable.pos(sel);
bsSur = 70;
motifMat=reshape(normProfile(acol(round(BStable.pos(sel))+[-bsSur:bsSur].*(BStable.dir(sel))),:), ...
    sum(sel),2*bsSur+1,width(normProfile));
[~,i] = sort(mean(motifMat,2),'descend');

%% show foot print
figure
subplot(10,1,1:2)
hold on
fill([40 101 101 40],[0 0 40 40],[254 198 1]./255)
fill([68 74 74 68],[0 0 40 40],[179 214 158]./255)
plot(movmean(mean(motifMat),5))
xticks(1:70:141)
xticklabels({'-70 bp','AGGGG','+70 bp'})

subplot(10,1,4:10)
imagesc(((motifMat(i,:))),[0 100])
colormap(brewermap(1000,'Blues'))
xticks(1:70:141)
xticklabels({'-70 bp','AGGGG','+70 bp'})
c=colorbar;
c.Location = 'southoutside'; 
%% show TSS
metaPro= metaProfilePromLenDivya(normProfile,'promEnd','position','afterTss',1000,'promLen',promoterLengthsORF);
[~,idx] = sort(sumProm,'descend','MissingPlacement','last');
metaPro1 = metaPro;
%%
figure
subplot(10,1,1:2)
hold on
fill([0 3033 3033 0],[0 0 2 2],[254 198 1]./255,'EdgeColor','none','FaceAlpha',0.5)

metaPro1(isnan(metaPro))=0;
plot(movmean(mean(metaPro1),100))
xticks([2333,3033,3733])
xticklabels({'-0.7kb','TSS','+0.7kb'})
xlim([2333 3733])
% ylim([0 2])
subplot(10,1,4:10)
metaPro(isnan(metaPro))=-1;
imagesc((metaPro(idx,:)))

colormap([.5 .5 .5;brewermap(1000,'Blues')])
xticks([2333,3033,3733])
xticklabels({'-0.7kb','TSS','+0.7kb'})
xlim([2333 3733])
caxis([-1 10])
ylim([1 5358])
c=colorbar;
c.Location = 'southoutside'; 
c.Label.String = 'ChEC signal';
%% B 
load('absForJing.mat')
load('BStables_1e3.mat')
load('promoterIDXvecFimp.mat')
BStable = BStables.MSN2;
for i = 1:width(nProfile)
    normProfile = nProfile(:,i);
    sumMat = zeros(4,1);
    bsSur = 30;
    sel = (BStable.type==1);

    motifMat=reshape(normProfile(acol(round(BStable.pos(sel))+[-bsSur:bsSur].*(BStable.dir(sel))),:),sum(sel),2*bsSur+1,width(normProfile));

    xTable.motifScore(i) = sum(motifMat(:));

    xTable.PromoterScore(i) = sum(normProfile(promoterIDXvecF==1),1);%AllpromoterMat
    xTable.mitoScore(i) = sum(normProfile(promoterIDXvecF==5),1);%mitochondira mat
    xTable.teloScore(i) = sum(normProfile(promoterIDXvecF==4),1);%TeloMat 
end
% clearvars -except xTable
%% B corr all measurements
figure
subplot(1,2,1)
scatter(log(xTable.PromoterScore),log(xTable.aro80Sig),30,'filled','MarkerFaceAlpha',0.5)
xlabel('Aro80 Spike-in')
ylabel('Motif measurement')
subplot(1,2,2)
scatter(log(xTable.motifScore),log(xTable.aro80Sig),30,'filled','MarkerFaceAlpha',0.5)
xlabel('Aro80 Spike-in')
ylabel('Motif measurement')

figure
allval= [xTable.PromoterScore,xTable.motifScore,xTable.mitoScore,xTable.aro80Sig];
labels = {'Promoter','Motifs','Aro80 spike-in','Mitochondria'};
[c,p]= corr(allval);
imagesc(corr(allval),[-1 1])
tempcmap = brewermap(100,'-RdBu');

colormap(tempcmap(26:75,:))
set(gca,'XTick',1:4,'YTick',1:4,'XTickLabel',labels,'YTickLabel',labels)
hold on
xline([0.5:1:3.5],'w')
yline([0.5:1:3.5],'w')
colorbar

%% C show corr. between all measurements
Allmeasure = log([datafiles_rp.AbsMotif; datafiles_rp.AbsPromoter;
    datafiles_rp.AbsMitochondria]');
validsel = find(all(Allmeasure,2)&all(Allmeasure>10,2)&all(Allmeasure<20,2));

Allnames = {'Motif','Promoter','Mitochondria'};
temppair = nchoosek(1:3,2);
clf
for i =1:height(temppair)
    subplot(1,3,i)
    scatter((Allmeasure(validsel,temppair(i,1))),(Allmeasure(validsel,temppair(i,2))),...
        10,'blue','filled','MarkerFaceAlpha',0.5)
    xlabel(Allnames{temppair(i,1)})
    ylabel(Allnames{temppair(i,2)})
    if i >1
        ylim([10 20])
    end
    [c,p] = corr(Allmeasure(validsel,temppair(i,1)),Allmeasure(validsel,temppair(i,2)));
    title(sprintf("r=%.2f, -log10(p)=%.2f",c,-log10(p)))
end

sgtitle('Compare all absolute values measurements')
%% D corr. abs. measurement in this study
figure
allval= [strains.meanAbs(strains.meanAbs(:,1)~=0&~isnan(strains.meanAbs(:,1)),[1:3])];
labels = {'Promoter','Motifs','Mitochondria'};
[c,p]= corr(allval);
imagesc(corr(allval),[-1 1])
tempcmap = brewermap(100,'-RdBu');

colormap(tempcmap(26:75,:))
set(gca,'XTick',1:4,'YTick',1:4,'XTickLabel',labels,'YTickLabel',labels)
hold on
xline([0.5:1:3.5],'w')
yline([0.5:1:3.5],'w')
colorbar
%% E diff. between repeats
% Absolute difference between any pair of measurements
for i = 1:3
    absdiffAny{i} = nonzeros(triu((log(Allmeasure(validsel,i)') - log(Allmeasure(validsel,i))), 1));
end

% Absolute difference between any pair of intra-repeat measurements
absdiffRp = cell(1,3);
for i = [1:699,height(strains)]
    if ~isnan(strains.meanAbs(i,1))
    tempvalues = ([strains.goodAbsSum{i};strains.badAbsSum{i}]);
    n=1;
    for measi = 1:3
        absdiffRp{measi} = [absdiffRp{measi};nonzeros(triu((log(tempvalues(:,measi)') - log(tempvalues(:,measi))), 1))];
    end
    end
end
%% E plot diff. between repeats distribution -- histogram
figure
hold on
histogram(abs(absdiffAny{2}),'Normalization','pdf','FaceColor',[166 206 226]./255,'EdgeColor','none')
histogram(abs(absdiffRp{2}),'Normalization','pdf','FaceColor',[54 160 71]./255,'EdgeColor','none')
legend('All promoter abs.','Intra-repeats')

xlabel('Log-scale absolute difference between any pair of measurements')
ylabel('Probability')
% set(gca,'xscale','log')
title(Allnames{i})


%% F distribution of good values and bad values
figure
rpsdrawlist = drawlist(1:25,:);
rpsdrawlist(maxcorr(1:25,:)<0.9) = 0;
gri = 13;
sel1 = nonzeros(acol(rpsdrawlist(:,[14:20])));
sel2 = nonzeros(rpsdrawlist(:,gri));
sel1 = strains.max_corr>=0.9;
sel2 = strains.max_corr<0.9;
Gv = cat(1,measureCollect{sel1});
Bv = cat(1,measureCollect{sel2});
i=2
hold on
h1=histogram(log(nonzeros(Gv(:,i))),'Normalization','pdf','FaceColor','#c49a6c','FaceAlpha',0.6,'EdgeColor','none');

h2=histogram(log(nonzeros(Bv(:,i))),'Normalization','pdf','FaceColor','#939598','FaceAlpha',0.6,'EdgeColor','none');
[~,b] = ttest2((nonzeros(Gv(:,i))),(nonzeros(Bv(:,i))));
if i==2
xlabel('log(Measurements)')
end
if i==1
ylabel('Probability')
end
title(sprintf('%s, -log10(p)= %.2f',Allnames{i},-log10(b)))
if i==3
legend('All other mutants',strrep(allgrammarlist{gri},'_',' '))
end

