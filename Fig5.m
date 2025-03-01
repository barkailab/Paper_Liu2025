load('PromoterGCcontent.mat')
%% select targets
change = {'Ldel','Fdel','FtoL','Fto2L'};
clearvars foldchangeMinMaxnF* foldchangeMinMaxnL corrBind* FLsel
for curj = 1:3
    for j = 1:2
    curInt = change{j};
    cursel = find(contains(denovoinfo.newtype,curInt)&denovoinfo.ref==curj);
    if ~isempty(cursel)
    % exclude bohdana strains
    cursel = cursel((strains.max_corr(denovoinfo.strainid(cursel)))>0.9);
    [~,idx] = sort(denovoinfo.denovoid(cursel));
%     if ~any(ismember(denovoinfo.newname(cursel),'Ref'))
%         cursel = [cursel; find(contains(denovoinfo.newname,['Ref',num2str(curj)]))];
%     end
     id0 = find(contains(denovoinfo.newname(cursel),{'F0','L0'}));
    if isempty(id0)
        cursel = [cursel; find(contains(denovoinfo.newname,[num2str(curj),'F0']))];
        id0 = find(contains(denovoinfo.newname(cursel),{'F0','L0'}));
    end
   corrF0 = corr(StrainSumProm(Protype<3,unique(denovoinfo.strainid(cursel(id0)))),StrainSumProm(Protype<3,denovoinfo.strainid(cursel)),'rows','pairwise');
    if j==1|j==2
    [~,refid] = min(corrF0);
    else
        refid = find(contains(denovoinfo.newname(cursel),['Ref',num2str(curj)]));
        refid = refid(1);
    end
    cur_sp = StrainSumProm(:,denovoinfo.strainid(cursel));
    log_sum = log2(cur_sp+700); % log for msn2 & my pertubations

    zscore= (log_sum-mean(log_sum,1,'omitnan'))./std(log_sum,[],1,'omitnan');
    z_log = zscore>3; 
    FLtargetList{curj,j} = find((sum(z_log,2)>=1)&(Protype<3)); 

    FLsel{curj,j} = cursel;
    if j ==1
        nLs = denovoinfo.nL(cursel);
        corrBindnLAll(Protype<3,curj) = corr(log_sum(Protype<3,:)',nLs,'rows','complete');
        foldchangeMinMaxnL(Protype<3,curj) = log_sum(Protype<3,refid)-log_sum(Protype<3,id0);
    elseif j ==2
        nFs = denovoinfo.nF(cursel);
        foldchangeMinMaxnF1(Protype<3,curj)= log_sum(Protype<3,refid)-log_sum(Protype<3,id0);
        foldchangeMinMaxnF2(Protype<3,curj) = log_sum(Protype<3,end)-log_sum(Protype<3,refid);
        corrBindnF1(Protype<3,curj) = corr(log_sum(Protype<3,id0:refid)',nFs(id0:refid),'rows','pairwise');
        corrBindnF2(Protype<3,curj) = corr(log_sum(Protype<3,refid:end)',nFs(refid:end),'rows','pairwise');
    elseif j ==3|j==4
        nFs = denovoinfo.nF(cursel);
        foldchangeMinMaxnFL(Protype<3,curj)= log_sum(Protype<3,refid)-log_sum(Protype<3,id0);
        corrBindnFL(Protype<3,curj) = corr(log_sum(Protype<3,:)',nFs,'rows','pairwise');

    end
    end
    end
end
%% wrtite to excel
for curj = 1:3
    RegList = table();
        cursel = cat(1,FLsel{curj,1:2});
    
    id0 = find(contains(denovoinfo.newname(cursel),{'F0','L0'}));
    corrF0 = corr(StrainSumProm(Protype<3,unique(denovoinfo.strainid(cursel(id0)))),StrainSumProm(Protype<3,denovoinfo.strainid(cursel)),'rows','pairwise');
    [~,refid] = min(corrF0);
    cur_sp = StrainSumProm(:,denovoinfo.strainid(cursel));
    log_sum = log2(cur_sp+700); % log for msn2 & my pertubations
    [~,orderU] = sort(log_sum(FLsubtar{curj,1},refid),'descend');
    [~,orderD] = sort(log_sum(FLsubtar{curj,2},refid),'descend');
    targetList = FLsubtar{curj,2}(orderD);
    %FLsubtar{curj,2}(orderD)];
    nHyd = denovoinfo.nL(cursel)+denovoinfo.nF(cursel);
    RegList.geneNames = GP.gene_infoR64.name(targetList);
    RegList.graduality_L = corrBindnLAll(targetList,curj);
    RegList.level_foldchange_L = foldchangeMinMaxnL(targetList,curj);
    RegList.graduality_F1 = corrBindnF1(targetList,curj);
    RegList.level_foldchange_F1 = foldchangeMinMaxnF1(targetList,curj);
    RegList.graduality_F2= corrBindnF2(targetList,curj);
    RegList.level_foldchange_F2 = foldchangeMinMaxnF2(targetList,curj);
    for j = 1:numel(nHyd)
        RegList.([denovoinfo.newname{cursel(j)},'_log2binding'])=log_sum(targetList,j);
    end
    writetable(RegList,'/home/labs/barkailab/jingliu/Documents/Finalfiles/Sfiles/FLGeneBinding.xlsx','Sheet',[num2str(curj),'down'])
end
%% plot Fig5B
%% scatter to show 
figure
for i = 1:3
    %seltargets
    seltar = unique(cat(1,FLtargetList{i,:}));
    subplot(3,3,3*(i-1)+1)
    hold on
    subsetseltarU = find(corrBindnF1(seltar,i)>0.75&corrBindnLAll(seltar,i)>0.75);
    subsetseltarD = find(corrBindnF1(seltar,i)<-0.5&corrBindnLAll(seltar,i)<-0.5);
    scatter(corrBindnF1(seltar,i),corrBindnLAll(seltar,i), 20,foldchangeMinMaxnL(seltar,i),'filled')
    [c,p]=corr(corrBindnF1(seltar,i),corrBindnLAll(seltar,i));
    [x,y] = polyfit(corrBindnF1(seltar,i),corrBindnLAll(seltar,i),1);

    x_fit = linspace(-1,1,100);
    y_fit = polyval(x, x_fit);
    plot(x_fit, y_fit, '-', 'color',[.5 .5 .5],'LineWidth', 2);
    xline([-0.5,0.75])
    yline([-0.5,0.75])

    FLsubtar{i,1} = seltar(subsetseltarU);
    FLsubtar{i,2} = seltar(subsetseltarD);
    colormap(brewermap(1000,'-RdBu'))           
    caxis([-3 3])
    xlabel('F left')
    ylabel('L')
    title(sprintf('color: fold change in L; r=%.2f, p=%.2f',c,-log10(p)))
    subplot(3,3,3*(i-1)+2)
    % get biggest change
    [~,id]=max([abs(foldchangeMinMaxnF1(seltar,i)),abs(foldchangeMinMaxnF2(seltar,i))],[],2);
    fdF = [foldchangeMinMaxnF1(seltar,i),foldchangeMinMaxnF2(seltar,i)];
    selfdF=[];
    for k = 1:height(fdF)
        selfdF(k) = fdF(k,id(k));
    end
    hold on
    scatter(corrBindnF2(seltar,i),corrBindnLAll(seltar,i), 20,selfdF,'filled')
    colormap(brewermap(1000,'-RdBu'))
       caxis([-3 3])

  
    xlabel('F right')
    ylabel('L')
    title('color: max fold change in F')
    subplot(3,3,3*(i-1)+3)
    hold on
    scatter(corrBindnF1(seltar,i),corrBindnF2(seltar,i), 20,selfdF,'filled')
    colormap(brewermap(1000,'-RdBu'))
        xlabel('F left')
    ylabel('F right')
    title('color: max fold change in F')
    caxis([-3 3])
    subsetseltarU = find(corrBindnF1(seltar,i)<-0.75&corrBindnF2(seltar,i)>0.5);
    subsetseltarD = find(corrBindnF1(seltar,i)>0.75&corrBindnF2(seltar,i)<-0.5);
    xline([-0.75,0.75])
    yline([-0.5,0.5])

    FLsubtar{i+3,2} = seltar(subsetseltarU);
    FLsubtar{i+3,1} = seltar(subsetseltarD);
end
%% plot Fig5C
change = {'Ldel','Fdel'};
logDBD = log2(TF_reference.MSN2_d02_604+700);
FLmap = [175 210 176;255 221 153]./255;
% clearvars foldchangeMinMaxnF* foldchangeMinMaxnL corrBind*
for curj = 1:3
%     for j = 1:2

    cursel = cat(1,FLsel{curj,1:2});
    
    id0 = find(contains(denovoinfo.newname(cursel),{'F0','L0'}));
    corrF0 = corr(StrainSumProm(Protype<3,unique(denovoinfo.strainid(cursel(id0)))),StrainSumProm(Protype<3,denovoinfo.strainid(cursel)),'rows','pairwise');
    [~,refid] = min(corrF0);
    cur_sp = StrainSumProm(:,denovoinfo.strainid(cursel));
    log_sum = log2(cur_sp+700); % log for msn2 & my pertubations
    [~,orderU] = sort(log_sum(FLsubtar{curj,1},refid),'descend');
    [~,orderD] = sort(log_sum(FLsubtar{curj,2},refid),'descend');
    targetList = [FLsubtar{curj,1}(orderU);FLsubtar{curj,2}(orderD)];
    nHyd = denovoinfo.nL(cursel)+denovoinfo.nF(cursel);
%
% figure
clf
axes('position',[0.05 0.05 0.23 0.85])
imagesc(log_sum(targetList,1:numel(FLsel{curj,1})),[12,18])
yline(numel(FLsubtar{curj,1})+0.5,'LineWidth',2)
xline([1:numel(FLsel{curj,1})]+0.5,'Color',[.75 .75 .75],'LineWidth',1)

hold on
geneids1 = GP.gene_infoR64.name(targetList);
% geneids1(topDBDinids1) = cellfun(@(x) sprintf('\\textbf{* %s}',x),geneids1(topDBDinids1),'UniformOutput',false);
set(gca,'XTick',1:numel(FLsel{curj,1}),'XTickLabel', denovoinfo.newname(FLsel{curj,1}),...
    'Ytick',1:numel(geneids1),'YTicklabel',geneids1,'XTickLabelRotation',90,'YDir', ...
    'reverse','ticklabelinterpreter','none','colormap',brewermap(1000,'Greens'),'FontSize',6)
% colorbar
title(References{curj})
axes('position',[0.32 0.05 0.23 0.85])
imagesc(log_sum(targetList,numel(FLsel{curj,1})+1:end),[12,18])
yline(numel(FLsubtar{curj,1})+0.5,'LineWidth',2)
xline([1:numel(FLsel{curj,1})]+0.5,'Color',[.75 .75 .75],'LineWidth',1)
hold on
geneids1 = GP.gene_infoR64.name(targetList);
% geneids1(topDBDinids1) = cellfun(@(x) sprintf('\\textbf{* %s}',x),geneids1(topDBDinids1),'UniformOutput',false);
set(gca,'XTick',1:numel(FLsel{curj,2}),'XTickLabel', denovoinfo.newname(FLsel{curj,2}),...
    'Ytick',[],'YTicklabel',geneids1,'XTickLabelRotation',90,'YDir', ...
    'reverse','ticklabelinterpreter','none','colormap',brewermap(1000,'Greens'),'fontsize',6)
% colorbar
title(References{curj})
% ----------------------------
% DBD
axes('position',[0.56 0.05 0.02 0.85])
imagesc(logDBD(targetList),[10,18])
set(gca,'colormap',brewermap(1000,'Greens'))
xticks(1)
xticklabels('DBD')
yticks([])

% OPN
axes('Position',[0.59,0.05 0.02 0.85])
imagesc(opn.opnScore(targetList),[-3 3])
set(gca,'colormap',brewermap(1000,'-RdBu'),'xtick',1,'XTickLabel','OPN','ytick',[])
% TATA
axes('Position',[0.62,0.05 0.02 0.85])
imagesc(simplifiedtata(targetList)',[0 9])
TATAmap = [.75,.75,.75;230/255,85/255,13/255];
TATAmap = [TATAmap;repmat([1,1,1],7,1)];
set(gca,'colormap',TATAmap,'xtick',1,'XTickLabel','TATA','ytick',[])
% GCcontent
axes('Position',[0.65,0.05 0.02 0.85])
imagesc(GCzscore(targetList),[-3 3])
set(gca,'colormap',brewermap(1000,'-RdGy'),'xtick',1,'XTickLabel','GC%','ytick',[])

axes('Position',[0.72,0.15 0.12 0.7])
hold on
mean_bind = mean(log_sum(FLsubtar{curj,1},:));
std_bind = std(log_sum(FLsubtar{curj,1},:));
for kn = 1:numel(nHyd)
scatter(nHyd(kn),mean_bind(kn),125,denovoinfo.nF(cursel(kn))==0,'filled','MarkerEdgeColor',[.75 .75 .75],'MarkerFaceAlpha',nHyd(kn)./max(nHyd))
end
% errorbar(nHyd,mean_bind,std_bind)

plot(nHyd(1:numel(FLsel{curj,1})),mean_bind(1:numel(FLsel{curj,1})),'Color','#FFDD99')
plot(nHyd(numel(FLsel{curj,1})+1:end),mean_bind(numel(FLsel{curj,1})+1:end),'Color','#B0D2B0')

colormap(gca,FLmap)

xlabel('nHyd')
% set(gca,'xtick',nHyd)
ylabel('log2 binding')
title(['top Pos & nL'])

axes('Position',[0.87,0.15 0.12 0.7])
hold on
mean_bind = mean(log_sum(FLsubtar{curj,2},:));
std_bind = std(log_sum(FLsubtar{curj,1},:));
for kn = 1:numel(nHyd)
scatter(nHyd(kn),mean_bind(kn),125,denovoinfo.nF(cursel(kn))==0,'filled','MarkerEdgeColor',[.75 .75 .75],'MarkerFaceAlpha',nHyd(kn)./max(nHyd))
end% errorbar(nHyd,mean_bind,std_bind)

plot(nHyd(1:numel(FLsel{curj,1})),mean_bind(1:numel(FLsel{curj,1})),'Color','#FFDD99')
plot(nHyd(numel(FLsel{curj,1})+1:end),mean_bind(numel(FLsel{curj,1})+1:end),'Color','#B0D2B0')
scatter(150,mean(logDBD(FLsubtar{curj,2})),50,'r','filled')
colormap(gca,FLmap)


xlabel('nHyd')
% set(gca,'xtick',nHyd)
ylabel('log2 binding')
title(['top Pos & nL'])
exportgraphics(gcf,[num2str(curj),'LF2.pdf'],'ContentType','vector')
%     end
end