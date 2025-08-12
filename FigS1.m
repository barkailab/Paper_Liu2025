cur_sp = StrainSumProm(:,FLid(1:25));

log_sum = log2(cur_sp+700); % log for msn2 & my pertubations

zscore= (log_sum-mean(log_sum,1,'omitnan'))./std(log_sum,[],1,'omitnan');
z_log = zscore>3;
%% targets together
targetList = find((sum(z_log,2)>=1)&(Protype<3)); 
z_logTar = z_log(targetList,:);
load('GP')
geneList = GP.gene_infoR64.name(targetList);

log_sPT =  log_sum(targetList,:); 
%% ordering targets by each ones' targets
nClu =17;
idxnow = kmeans(log_sPT,nClu,'Distance','correlation','Replicates',100,'MaxIter',1000,'Options',statset('UseParallel',1,'UseSubstreams',0));
% order clusters
idxnow(idxnow==10)=12;
for i = 1:17
    temp(i)  =sum(idxnow==i);
end
[~,orderidxnow]=sort(temp,'descend');

ids = idxnow;
for i = 1:max(idxnow)
    ids(idxnow==orderidxnow(i))=i;
end
idxnow=ids;
median_clus = arrayfun(@(x) median(log_sPT(idxnow==x,:)),1:max(idxnow),'UniformOutput',false);
median_clus = cat(1,median_clus{:});
corrpm = corr(median_clus',log_sPT');
corrTopCluster=corrpm'.*([1:1:nClu]==idxnow);

[tempTable,newIdx]=sortrows(table(idxnow,corrTopCluster));

figure

imagesc(corr(log_sPT(newIdx,:)'),[0,1])
colorbar
title(sprintf('correlate abs promoters(want %d clusters)', nClu))

%% FigS1 A show TFs
neworder = [];
counts = [];
for i = [12,9,2,7,3,5,16,15,13,6,1,8,17,14,4]
     [~,order] = sort(median(log_sPT(idxnow==(i),:),2),'descend'); 
     a = find(idxnow==(i));
     neworder = [neworder;a(order)];
     counts = [counts,numel(a)];
end
counts = cumsum(counts);
orderTar=neworder;
% order TFs
orderTF = [1 22 21 4 25 15 7 18 12 13 14 20 10 17 6 8 3 23 11 9 5 24 2 16 19];

figure
imagesc(log_sPT(orderTar,orderTF)',[12 17])

xline((counts)+0.5,'color','k','linewidth',1)

set(gca,'colormap',brewermap(1000,'Greens'),'Ytick',1:25,'YTickLabel',TF(orderTF'))
genenames = table();
genenames.Genename=GP.gene_infoR64.name(orderTar);
writetable(genenames,'TableS5_SuppValues.xlsx','Sheet','FigS1A')
%% FigS1B
%% conservations
% these are just MUSCLE aligned results
load('int_FL_ortho.mat')
load('info_ortho.mat')
%% conservation by blosum
blomat = blosum62;
blomat(25,:)=-4;
blomat(:,25)=-4;

for i = 1:height(TFinfo)-1
    if i >=21
        c = i+1;
    else 
        c=i;
    end 
    intmatrix = FLintmat{c};
    scerseq = othoinfo{c}.type ==1;
    indmatrix = sub2ind(size(blomat),intmatrix,repmat(intmatrix(scerseq,:),height(intmatrix),1));
    scorematrix = blomat(indmatrix);
    introw = contains(othoinfo{c}.Var1,{'skdz','cgla','kafr','tbla','zrou','tglo','klac','kwal'});
    avgvector = mean(scorematrix(introw,:),1);
    scerconv = avgvector(intmatrix(scerseq,:)~=25);

    TFinfo.scerConv{i} =scerconv;
end

%% get nonDBD conv and DBD conv
for i = 1:25
    seqfile = fastaread(['/home/labs/barkailab/jingliu/Documents/Finalfiles/orthologs/',TFinfo.TF{i},'_orthologs.fasta']);
    seqFL = seqfile(contains({seqfile.Header},'scer_A')).Sequence;

   if TFinfo.nonDBD_start(i)<10
        nonDBDst = 1;
    else
        nonDBDst = TFinfo.nonDBD_start(i);
    end
    nonDBD = TFinfo.scerConv{i}([nonDBDst:TFinfo.nonDBD_end(i)]);

    if TFinfo.FL_rm(i)>0
        DBDsel = true(1,length(seqFL));        
        DBDsel(sort([nonDBDst:TFinfo.nonDBD_end(i),TFinfo.FL_rm(i):TFinfo.FL_rmend(i)]))=false;
        DBD =  TFinfo.scerConv{i}(DBDsel);
        FLsel = true(1,length(seqFL));        
        FLsel(TFinfo.FL_rm(i):TFinfo.FL_rmend(i))=false;
      
        trueFL =  TFinfo.scerConv{i}(FLsel);
    else
        DBDsel = true(1,length(seqFL));        
        DBDsel(sort(nonDBDst:TFinfo.nonDBD_end(i)))=false;
        DBD =  TFinfo.scerConv{i}(DBDsel);
        trueFL  =  TFinfo.scerConv{i};
    end
    TFinfo.nonDBDConv(i) = mean(nonDBD);
    TFinfo.DBDconv(i) = mean(DBD);
end
%% FigS1B scatter of conservation
figure
scatter(TFinfo.nonDBDConv,TFinfo.DBDconv,'filled')
xlim([-2.5,4])
ylim(xlim)
hold on
plot(xlim,xlim)
xlabel('NonDBD')
ylabel('Respective DBD')
title('Conservation')

%% Fig S1I
