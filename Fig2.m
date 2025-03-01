
%% load oligos
nonDBDs = fastaread('nonDBD_seq.fasta');
load('FirstBatch.mat')
%% give AA names
for i = 1:height(FirstBatch)-1
frq = accumarray(FirstBatch.aaint{i}',1,[20,1]);
name = 'MSN2';
scaffold = [];
hydaa = [];
    for j = aaorder
        j
        if frq(j)>0
           name = [name,'_a',int2aa(j),num2str(frq(j))]
        end
        if j==aaorder(7)
            scaffold = name;
        end
        if find(ismember(aaorder(end-5:end),j))
            if frq(j)>0
                hydaa = [int2aa(j),num2str(frq(j))];
            end
        end
    end

    FirstBatch.newname{i} = name;
    FirstBatch.scaffold{i} = scaffold(5:end);
    FirstBatch.hydaa{i} = hydaa;
    FirstBatch.frq{i} = frq;
end
FirstBatch.frq{end} = zeros(20,1);
FirstBatch.newname(27:28) = {'MSN2','MSN2_d02_604'};
%% getting strainsid 
for i = 1:height(FirstBatch)
    strainid = min(find(ismember(strains.strain,FirstBatch.newname{i})));  
    if isempty(strainid)
         strainid = find(ismember(strains.strain,FirstBatch.newname{i}));
    end 
    FirstBatch.maxcorr(i)= strains.max_corr(strainid);
    FirstBatch.strainid(i) = strainid;
end
%% corrWT corrDBD
MDcorrWT = corr(StrainSumProm(Protype<3,FIrstBatch.strainid(FIrstBatch.Order>0)),TF_reference.MSN2(Protype<3),'rows','pairwise');
MDcorrDBD = corr(StrainSumProm(Protype<3,FIrstBatch.strainid(FIrstBatch.Order>0)),TF_reference.MSN2_d02_604(Protype<3),'rows','pairwise');
%% Fig2A bottom plot sequence
figure
tempsel = 22;
AAseq =FIrstBatch.aaint{tempsel};
imagesc(AAseq,[1,20])
set(gca, 'ytick',[],'xtick',[])
colormap(AAcolormap)
yline(1.5,'k')
ylim([0.5,1.5])
xlim([0,numel(FirstBatch.aaint{tempsel})])
title(strrep(FIrstBatch.newname{i},'_a',' '))
AAcompmat = cat(2,FIrstBatch.frq{:});
AAcompmat(21,:) = AAcompmat(4,:)+AAcompmat(7,:);
AAcompmat(22,:) = AAcompmat(2,:)+AAcompmat(12,:);
AAcompmatO = AAcompmat(:,[FIrstBatch.ordernew(1:end-2);27]);
%%  get colormap for each AA
cmap =[];

for i = 1:20
    if ismember(i,[3 6 16 17])

% Define the starting and ending colors
endColor =[227 153 214]/255;
    else
        endColor = AAcolormap(i,:);
    end
if ismember(i,aa2int('FY'))
    endColor = AAcolormap(i,:)*83/98;
elseif i == aa2int('L')
    endColor = [224 197 2]/255;
end

startColor = [1, 1, 1];  % White

% Define the number of colors in the colormap
numColors = max(AAcompmatO(i,:)); 

% Create a colormap matrix with gradual transition
cmap{i} = [linspace(startColor(1), endColor(1), numColors)', ...
        linspace(startColor(2), endColor(2), numColors)', ...
        linspace(startColor(3), endColor(3), numColors)'];
end
%% add light one for DE
AAcompmatO(21,:) = AAcompmatO(4,:)+AAcompmatO(7,:);


endColor =[0.97 0.85 0.9]

startColor = [1, 1, 1];  % White

% Define the number of colors in the colormap
numColors = 100;  % Adjust as needed for smoothness

% Create a colormap matrix with gradual transition
cmap{21} = [linspace(startColor(1), endColor(1), numColors)', ...
        linspace(startColor(2), endColor(2), numColors)', ...
        linspace(startColor(3), endColor(3), numColors)'];
%% add light one for KR
AAcompmatO(22,:) = AAcompmatO(12,:)+AAcompmatO(2,:);

endColor =[0.68 0.85 0.99];

startColor = [1, 1, 1];  % White

% Define the number of colors in the colormap
numColors = 20;  % Adjust as needed for smoothness

% Create a colormap matrix with gradual transition
cmap{22} = [linspace(startColor(1), endColor(1), numColors)', ...
        linspace(startColor(2), endColor(2), numColors)', ...
        linspace(startColor(3), endColor(3), numColors)'];
%% Order scaffolds
[tempscaffolduni,~,tempid] = unique(FIrstBatch.scaffold(1:26));
tempnewid = tempid;
tempnewid (tempid==5) = 1;
tempnewid (tempid==4) = 2;
tempnewid (tempid==6) = 3;
tempnewid (tempid==9) = 4;
tempnewid (tempid==1) = 5;
tempnewid (tempid==7) = 6;
tempnewid (tempid==2) = 7;
tempnewid (tempid==3) = 8;
tempnewid (tempid==8) = 9;
FIrstBatch.scaffoldidnew(1:26) = tempnewid;
order = [];
for i = 1:9
    tempid = find(tempnewid ==i);
    [~,tempidorder] = sort(AAcompmat(21,tempid),'descend');
    order = [order;tempid(tempidorder)];
end
FIrstBatch.ordernew(1:26) = order;
%% Fig2B show my design
xticknames = {'N','Q','S','T','G','P','DE','KR','Hyd','Y','L'};

temp = brewermap(9,'Paired');
scaffoldcp = temp;
scaffoldcp(9,:) = temp(7,:);
scaffoldcp(7,:) = temp(9,:);
scaffoldcp(3,:) = temp(4,:);
scaffoldcp(4,:) = temp(3,:);
figure

c=0;

selid = 0;
for i = [aa2int('NQSTGP'),21,22,aa2int(['FYL'])]
     c=c+1;
     temph = 1;
  
     for j = 1:9
         tempid = find(FIrstBatch.scaffoldidnew==j);
    [~,tempidorder] = sort(AAcompmat(21,tempid),'descend');
    selid = tempid(tempidorder);
         temph = temph - 0.02*numel(tempid)-0.015;  
         if i ==3
         axes('Position',[0.06,temph,0.03,0.02*numel(tempid)])
          imagesc(ones(1,numel(tempid)),1:numel(tempid),repmat(j,1,numel(tempid)),[1 9])
             set(gca, 'colormap',scaffoldcp,'xtick',[],'ytick',[])

         end
 
    if ismember(i,[21,22,aa2int('FYL')])
        axes('Position',[0.05*(c+1)+0.02*(c-6),temph,0.05,0.02*numel(tempid)])
    else
         axes('Position',[0.05+0.05*c,temph,0.05,0.02*numel(tempid)])
    end
    imagesc(ones(1,numel(tempid)),1:numel(tempid),AAcompmat(i,selid)',[0,max(AAcompmat(i,:))])
    set(gca, 'colormap',cmap{i},'XTick',[],'ytick',[])
  
     nums =string(AAcompmat(i,selid)');
    nums(strcmp(nums,"0"))=" ";
    if ismember(i,[21,22,aa2int('FYL')])
        withchange = find(AAcompmat(i,selid)~=0);
        for j = withchange
            rectangle('Position',[0.5 j-0.5 1 1],'EdgeColor','k','FaceColor','none')
        end
    else
    yline([0.5:1:numel(tempid)-0.5],'k')
     end
  
    if i==aa2int('F')
        nums = FIrstBatch.hydaa(selid);
       
    end
    text(repmat(0.6,1,numel(tempid)),1:numel(tempid),nums)
     end
      set(gca, 'colormap',cmap{i},'XTick',1,'XTickLabel',xticknames{c},'ytick',[],'xticklabelrotation',0)
  
end
exportgraphics(gcf,'Fig2B showdesign.eps','ContentType','vector')
%% Fig 2C
figure
hold on
% no hydrophobic no DE
sel = ((AAcompmat(21,1:26)==0)' & cellfun(@isempty,FIrstBatch.hydaa(1:26)));
scatter(MDcorrDBD(sel),FIrstBatch.maxcorr(sel),50,FIrstBatch.scaffoldidnew(sel),'filled','marker','square','MarkerEdgeColor','k')
% no hydrophobic and yesDE
sel = ((AAcompmat(21,1:26)~=0)' & cellfun(@isempty,FIrstBatch.hydaa(1:26)));
scatter(MDcorrDBD(sel),FIrstBatch.maxcorr(sel),50,FIrstBatch.scaffoldidnew(sel),'filled','marker','diamond','MarkerEdgeColor','k')
%  hydrophobic and noDE
sel = ((AAcompmat(21,1:26)==0)' & ~cellfun(@isempty,FIrstBatch.hydaa(1:26)));
scatter(MDcorrDBD(sel),FIrstBatch.maxcorr(sel),50,FIrstBatch.scaffoldidnew(sel),'filled','marker','^','MarkerEdgeColor','k')

%  hydrophobic and yesDE
sel = ((AAcompmat(21,1:26)~=0)' & ~cellfun(@isempty,FIrstBatch.hydaa(1:26)));
scatter(MDcorrDBD(sel),FIrstBatch.maxcorr(sel),50,FIrstBatch.scaffoldidnew(sel),'filled','marker','o','MarkerEdgeColor','k')

    colormap(gca,scaffoldcp)
    yline(0.9,'color',[.5 .5 .5],'LineStyle','--')
    text(0,0.87,'Reproducible')
hold on
% scatter(corrDBD(15,1),maxcorr(15,1),70,'^','filled')
scatter(0.582,0.9966,70,'r','pentagram','filled')
legend({'Scaffold only','Scaffold+DE','Scaffold+Hyd','Scaffold+DE+Hyd'})
xlabel('corrDBD')
ylabel('Reproducibility')

%% Order reproducible denovo IDR - Msn2DBD
orderClust = [6 5 19 3 13 8 7 11 12 17 1 14 15 9 10 16 2 20 21 18 4];
goodids = [FIrstBatch.strainid(FIrstBatch.maxcorr>0.9)];
names = [repmat({' '},1,24),'Msn2','DBD'];
goodsel =find(FIrstBatch.maxcorr>0.9);
%% calculate in vitromotif binidng
BStable = BStables.(TF{15}); % 15 is Msn2
bsSur = 150;
sel = (BStable.type==1);
clearvars motifMat
normProfile = [];
c=1;
for j = goodids(orderClust)'
    tempid = find(ismember({DK_normP.name}, [strains.strain{j},'.mat']));
    if isempty(tempid)
        oldnameid = find(FIrstBatch.strainid==j);
        tempid = find(ismember({DK_normP.name}, [FIrstBatch.Var1{oldnameid},'.mat']));
    end
    if numel(tempid)>1
        tempid = tempid(1);
    end
    temp = load([DK_normP(tempid).folder '/' DK_normP(tempid).name]); 
    if isfield(temp, 'normProfile')
    normProfile(:,c) = temp.normProfile;
    elseif isfield(temp,'medianNorm')
    normProfile(:,c) = temp.medianNorm;    
    end
    c=c+1;
end
motifMat=reshape(normProfile(acol(round(BStable.pos(sel))+[-bsSur:bsSur].*(BStable.dir(sel))),:),sum(sel),2*bsSur+1,width(normProfile));
mean_promoter = squeeze(mean(motifMat(:,35:105,:),[1,2])); % +- 35 bp regions
promBaseBind=mean(normProfile(promoterIDXvecF==1,:),1);
bindScore=mean_promoter./promBaseBind';
strains.bindScore(goodids(ordernow)) = bindScore;
%% Fig2D
%% heatmap - corr. matrix
ordernow = orderClust;
figure
subplot(1,20,2:18)
imagesc(corr(StrainSumProm(Protype<3, goodids(ordernow)),'rows','pairwise'),[0,1])
set(gca,'ytick',1:numel(ordernow),'yticklabel',names(ordernow),'xtick',[])
colormap(gca,brighten(flipud(brewermap(1000,'-Blues')),0.3))
hold on
xline([0.5:1:numel(ordernow)-0.5],'w','linewidth',0.1)
yline([0.5:1:numel(ordernow)-0.5],'w','linewidth',0.1)
caxis([0.5816,1])
subplot(1,20,1)
imagesc(FIrstBatch.scaffoldidnew(goodsel(ordernow)),[0,9])
set(gca,'colormap',[1 1 1;scaffoldcp],'YTick',[])
axis equal
xlim([0.5,1.5])  
xticks([])
subplot(1,20,19)
% 'M' get Motif
imagesc(log2(strains.bindScore(goodids(ordernow))),[0 4])
colormap(gca,brewermap(1000,'Greens'))
yticks([])
colorbar
xticks([])
subplot(1,20,20)
% 'O' get OPN
cur_sp = StrainSumProm(:,goodids(ordernow));
log_sum = log2(cur_sp+700); % supress noise
zscore= (log_sum-mean(log_sum,1,'omitnan'))./std(log_sum,[],1,'omitnan');
z_log = zscore>3;
for i = 1:width(z_log)
    tempopn(i) = median(opn.opnScore(find(z_log(:,i))),'omitnan');
end
imagesc(tempopn',[-1 1])
colormap(gca,brewermap(1000,'-RdBu'))
yticks([])
colorbar
xticks([])
%% top and bottom bar corrMsn2, corr Msn2DBD
clf
subplot(20,1,1)
imagesc(corr(StrainSumProm(Protype<3, goodids([22])),StrainSumProm(Protype<3, goodids(ordernow)),'rows','pairwise'),[0,1])
set(gca,'ytick',1:numel(ordernow),'yticklabel',names(ordernow),'xtick',[])
colormap(gca,brighten(flipud(brewermap(1000,'-Blues')),0.3))
hold on
xline([0.5:1:numel(ordernow)-0.5],'w','linewidth',0.1)
% yline([0.5:1:numel(ordernow)-0.5],'w','linewidth',0.1)
thresh = round(corrWT(15,2),1);
clim([thresh,1])
subplot(20,1,3)
imagesc(corr(StrainSumProm(Protype<3, goodids([23])),StrainSumProm(Protype<3, goodids(ordernow)),'rows','pairwise'),[0,1])
set(gca,'ytick',1:numel(ordernow),'yticklabel',names(ordernow),'xtick',[])
colormap(gca,brighten(flipud(brewermap(1000,'-Blues')),0.3))
hold on
xline([0.5:1:numel(ordernow)-0.5],'w','linewidth',0.1)
% yline([0.5:1:numel(ordernow)-0.5],'w','linewidth',0.1)
clim([thresh,1])
%% Left: plot AA composition
figure
% Apply the colormap
xticknames = {'N','Q','S','T','G','P','DE','KR','Hyd','Y','L'};

c=0;
for i = [aa2int('NQSTGP'),21,22,aa2int(['FYL'])]
    c=c+1;
    axes('Position',[0.05+0.05*c,0.1,0.05,0.8])
    imagesc(ones(1,numel(ordernow)),1:numel(ordernow),goodComMat(i,ordernow)',[0,max(goodComMat(i,:))])
    set(gca, 'colormap',cmap{i},'XTick',1,'XTickLabel',xticknames{c},'ytick',[])
    nums =string(goodComMat(i,ordernow)');
    nums(strcmp(nums,"0"))=" ";
    if i==aa2int('F')
        nums = FIrstBatch.hydaa(goodsel(ordernow));
       
    end
    text(repmat(0.6,1,numel(ordernow)),1:numel(ordernow),nums)
    yline([0.5:1:numel(ordernow)-0.5],'k')
    ylim([0.5 numel(ordernow)+0.5])
end
exportgraphics(gcf,'Fig2DcorrBind.eps','ContentType','vector')
%% Figure S2B
