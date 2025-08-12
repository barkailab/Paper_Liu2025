% get targets for each
allFullReplaced = cellfun(@(x) x(end),pertubationinfo.cursels(1:21));
selTargets = arrayfun(@(x) find(z_log(:,x)), allFullReplaced,'UniformOutput',0);
perbOrder = bestorder(cur_sp(Protype<3, allFullReplaced)');
refcolormap = [166 206 226;1 1 1;246 153 153]./255;
%% include felix strains
AllTF = load('allTFs.mat');
DivyaDBD = load('medianSumPromNewAll.mat');
descTF = readtable('descTF.xlsx');
toadd = find(~ismember(TFdescription.TF,descTF.TF));
tempH = height(descTF);
n=0;
for i = toadd'
    n=n+1;
    descTF.TF{tempH+n} = TFdescription.TF{i};
    descTF.DBDname{tempH+n} = TFdescription.DBD{i};
end
%% get all DBDs
DBDsumProm = [];
listDBD = descTF.DBDname(find(cellfun(@(x) ~isempty(x),descTF.DBDname)));
for i = 1:numel(listDBD)
    try
        DBDsumProm(:,i) = DivyaDBD.medianSumPromNewAll.(listDBD{i});
    catch
        DBDsumProm(:,i) = StrainSumProm(:,find(ismember(strains.strain,listDBD{i})));
    end
end
% DBDsumProm(:,47) = DivyaDBD.medianSumPromNewAll.(TFdescription.dkDBD{25});
%% get all TFs
toadd = find(~ismember(TFinfo.TF,AllTF.currStrainsFelix));
listTF =  [AllTF.currStrainsFelix(find(~contains(AllTF.currStrainsFelix,'DBD')));TFinfo.TF(toadd)];

TFsumProm = [];
for i =196:numel(listTF)
    try 
    TFsumProm(:,i)= AllTF.sumPromFelix(:,find(ismember(AllTF.currStrainsFelix,listTF{i})));
    catch
        if i ==182|i==192
            TFsumProm(:,i) = DivyaDBD.medianSumPromNewAll.GIS1;
        else
        TFsumProm(:,i) =  StrainSumProm(:,find(ismember(strains.strain,listTF{i})));
        end
    end
end
%% get hugh and hahn
nonemty = find(cellfun(@(x) ~isempty(x),HahnFeatures.coactivator));
allHahnid = zeros(size(Protype));
[allHahn,~,allHahnid(nonemty)] = unique(HahnFeatures.coactivator(nonemty));
nonemtya= nonemty;
nonemty = find(cellfun(@(x) ~isempty(x),PughFeatures.FeatureClassLevel1));
[allPugh,~,allPughid(nonemty)] = unique(PughFeatures.FeatureClassLevel1(nonemty));
groups = [1,2,3,3,4,4,4,4,4];
groupedallPughid = zeros(size(Protype));
groupedallPughid(nonemty) = groups(allPughid(nonemty));
PughFeat = {'RP','STM','noTF','others'};

%% OPN
[~,opnorder] = sort(opn.opnScore(~isnan(opn.opnScore)),'descend');
temp = find(~isnan(opn.opnScore));
opnorder = temp(opnorder);
opnorder = opnorder(ismember(opnorder,nonemty));
opnorder = opnorder(ismember(opnorder,nonemtya));
figure
subplot(1,20,1)
imagesc(opn.opnScore(opnorder),[-3 3])
colormap(gca,brewermap(1000,'-RdBu'))
xticks([])
subplot(1,20,3:7)
zallDBD =  (DBDsumProm-mean(DBDsumProm,1,'omitnan'))./std(DBDsumProm,[],1,'omitnan');
imagesc(zallDBD(opnorder,:),[0 3])
colormap(gca,brewermap(1000,'Greens'))
yticks([])
xticks([])
subplot(1,20,8:12)
zallTF =  (TFsumProm-mean(TFsumProm,1,'omitnan'))./std(TFsumProm,[],1,'omitnan');
imagesc(zallTF(opnorder,:),[0 3])
colormap(gca,brewermap(1000,'Greens'))
yticks([])
xticks([])
clearvars zallDBD 
subplot(1,20,13:17)
order = bestorder(cur_sp(Protype<3, allFullReplaced)');
imagesc(zscore(opnorder,cat(1,pertubationinfo.cursels{order})),[0 3])
colormap(gca,brewermap(1000,'Greens'))
yticks([])
xticks([])
subplot(1,20,19)
imagesc(groupedallPughid(opnorder),[0,4])
colormap(gca,[230 230 230;211 157 85;141 11 65 ; 214 207 180; 255 248 230]./255)
yticks([])
subplot(1,20,20)
imagesc(allHahnid(opnorder),[0,2])
colormap(gca,[230 230 230;39 68 93;113 187 178]./255)
yticks([])
%% get sets
Msn2set = find(zDBD(:,6)>3&Protype<3);
DBDset = find(zDBD(:,5)>3&Protype<3);
for i = 1:21
    PerbTarSet{i}=find(any(zscore(:,pertubationinfo.cursels{i})>3,2));
end 
%% Fig7B scatter of targets and OPN:

% TATA
ar= simplifiedtata;
mTATA=sum(ar(Msn2set)~=0)./numel(Msn2set);
DBDTATA=sum(ar(DBDset)~=0)./numel(Msn2set);

for i = 1:21
    ratioTATA(i) = sum(ar(PerbTarSet{i})~=0)./numel(PerbTarSet{i});
end
figure
hold on
% histogram(ratioTATA)
xline(mTATA)
xline(DBDTATA)
[~,order] = sort(ratioTATA,'ascend');
scatter(ratioTATA(order),1:21,30,pertubationinfo.ref(1:21),'filled')
set(gca,'ytick',1:21,'yTickLabel',pertubationinfo.fulltype(order),'ydir','reverse')
xlabel('median TATA%')
ar= GCzscore;% or OPN
mTATA=median(ar(Msn2set));
DBDTATA=meidan(ar(DBDset));
for i = 1:21
    temp(i) = median(ar(PerbTarSet{i}));
end
figure
hold on
xline(mTATA)
xline(DBDTATA)
[~,order] = sort(temp,'ascend');
scatter(temp(order),1:21,30,pertubationinfo.ref(1:21),'filled')
colormap(rfcolormap)
set(gca,'ytick',1:21,'yTickLabel',pertubationinfo.fulltype(order),'ydir','reverse')
xlabel('median TATA%')


%% Fig7D correlate full replacements
orderP = bestorder(cur_sp(Protype<3, allFullReplaced)');
figure
imagesc(corr(cur_sp(Protype<3,allFullReplaced(orderP))),[0 1])
set(gca,'ytick',1:21','YTickLabel',pertubationinfo.fulltype(order),'xtick',1:21','xTickLabel',pertubationinfo.fulltype(order))
colormap(gca,brewermap(1000,'Blues'))
caxis([0.5 1])
%% Fig7D bottom
figure
distinctA = 6;% these values for which one you wanna look at
distinctB = 10;
targets = unique(cat(1,PerbTarSet{[distinctA,distinctB]}));
scatter(log_sum(targets, allFullReplaced(distinctA)),log_sum(targets, allFullReplaced(distinctB)),30,'b','filled','MarkerFaceAlpha',0.5)
xlabel(pertubationinfo.fulltype{distinctA})
ylabel(pertubationinfo.fulltype{distinctB})
hold on
plot(xlim,xlim)

title(num2str(corr(log_sum(targets, allFullReplaced((distinctA))),log_sum(targets, allFullReplaced((distinctB))))))



%% intra-perb

for i = 1:21
    xn= denovoinfo.(pertubationinfo.var{i})(refsel(pertubationinfo.ref(i)))-denovoinfo.(pertubationinfo.var{i})(sel(pertubationinfo.cursels{i}));
    TempcorrRefs = corr(StrainSumProm(Protype<3,denovoinfo.strainid(refsel(pertubationinfo.ref(i)))),cur_sp(Protype<3,pertubationinfo.cursels{i}));
    pertubationinfo.all_xRef{i} = TempcorrRefs;
    pertubationinfo.worst_xRef(i) = min(TempcorrRefs);
    pertubationinfo.xn{i} = xn;
end
temp=1:21;
Compared = [];
marksize = [];
for i=1:21
    toCompare = cat(1,pertubationinfo.cursels{temp(temp~=i)});
    c=0;

    j=allFullReplaced(i);
    c=c+1;
    marksize{i,c} = repmat(pertubationinfo.xn{i}(c),1,numel(toCompare));
    Compared{i} = corr(cur_sp(Protype<3,j),cur_sp(Protype<3,toCompare));

end

%% Fig7E all dots
catCompare = Compared;
[x,order] = sort(cellfun(@median, catCompare));
order = orderP;
% figure
clf
hold on
for i = 1:21
    scatter(catCompare{order(i)},repmat(i,size(catCompare{order(i)}))+0.2*rand(1,numel(catCompare{order(i)})),30,refcolormap(pertubationinfo.ref(order(i)),:),'filled','MarkerFaceAlpha',0.3)
end
set(gca,'ytick',1:21,'yticklabel',pertubationinfo.fulltype(order),'yDir','reverse')

