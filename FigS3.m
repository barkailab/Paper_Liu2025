% fig.S3A
%% make normprofile of interested strains
clearvars motifMat
normProfile = zeros(12157105,1);
i=27; % choose TF id you are interested in
for j = 31:36  
    clear motifMat
    strainid = drawlist(i,j);
    if strainid>0
        if strains.bestrp(strainid)==0
            mid = find(ismember({DK_normP.name},[strains.strain{strainid},'.mat']));
            temp = load([DK_normP(mid).folder '/' DK_normP(mid).name]);
            if isfield(temp, 'medianNorm')
                normProfile(1:size(temp.normProfile,1),j-30) = temp.medianNorm;
            elseif isfield(temp, 'normProfile')
                normProfile(1:size(temp.normProfile,1),j-30) = temp.normProfile;
            end
        else
            temp = load([datafiles_rp(strains.bestrp(strainid)).folder '/' datafiles_rp(strains.bestrp(strainid)).name]);
            normProfile(1:size(temp.normProfile,1),j-30) = temp.normProfile;
        end
    end
end
%% map it to desired invitro motifs
    
figure
TFn=0;
for i= denovoTFid(1:9)'% here you choose the interested TFid for interested DBD-recoginized motifs 
    TFn=TFn+1;% fig number
    BStable = BStables.(TF{i});
    bsSur = 70;
    sel = (BStable.type==1);
    motifMat=reshape(normProfile(acol(round(BStable.pos(sel))+[-bsSur:bsSur].*(BStable.dir(sel))),:),sum(sel),2*bsSur+1,width(normProfile));
    XDwt{TFn} = motifMat;
    % get motif enrichment 
    motifBind=squeeze(mean(motifMat(:,35:105,:),[1,2]));
    promBaseBind=mean(normProfile(promoterIDXvecF==1,:));
    XDbindScore(TFn,1:numel(promBaseBind))=motifBind'./promBaseBind;
    motifMat = cat(3,DDwt{TFn}(:,:,1:2),motifMat);
    mean_motifMat = squeeze(mean(motifMat,1));
           
    for j = 1:6
        subplot(10,6,6*(TFn-1)+j)
        hold on
        plot(movmean(mean_motifMat(:,1),5),'--','linewidth',0.5,'Color','#E7AA77')
        plot(movmean(mean_motifMat(:,2),5),':','linewidth',0.5,'Color','#B8C4BB')
        
         plot(movmean(mean_motifMat(:,j+2),5),'m-','linewidth',0.75,'Color','#E9F5ED')
        xticks([])
        set(gca,'fontsize',4)
        xlim([1,141])
        %                 ylim([0 10])
        if find(i==orderTF)==1
            title(do_labels{j},'FontSize',12)
        end
        if j==1
            ylabel(TF{i},'Rotation',0,'FontSize',8)
        end
            xlimv = xlim;
        ylimv = ylim;
        bg = patch([xlimv(1) xlimv(2) xlimv(2) xlimv(1)],[ylimv(1) ylimv(1) ylimv(2) ylimv(2)], ...
        corrmap(round(maxcorr(27,j+30)*1000)-1,:),'FaceAlpha',1,'EdgeColor','none');
        uistack(bg,'bottom')

    end
end
exportgraphics(gcf,'invitromotif.eps','ContentType','vector')

%% fig. S3B seqlogo -FL
figure
TFn=0;
for i = denovoTFid(1:9)'
    TFn=TFn+1;
    subplot(1,9,TFn)
    temp = load([DK_normP(min(find(ismember({DK_normP.name},[TF{i} '.mat'])))).folder '/' DK_normP(min(find(ismember({DK_normP.name},[TF{i} '.mat'])))).name]);
    if isfield(temp, 'medianNorm')
        normProfile(1:size(temp.medianNorm,1),1) = temp.medianNorm;
    else
        normProfile(1:max(size(temp.normProfile)),1) = temp.normProfile;
    end
 
    [out,sc4red]=mer_occupancy(normProfile(:,1),7,'window',21);
%         in out you have 8192 7mers and binding scores
    [~,top10]=maxk(out.score,10);

    [finalFreq]= seqLogoMin(out.mers.seq(top10),out.score(top10),'showFigs',true);
    title(TF{i});
end 


%% denovo IDR-DBD fusions-logo
for i = denovoTFid(1:6)'
    figure
    for j = 31:36
        strainid = drawlist(i,j);
        
        if strains.bestrp(strainid)==0
            mid = min(find(ismember({DK_normP.name},[strains.strain{strainid},'.mat'])));
            if i==15
                try 
                    temp = load([DK_normP(mid).folder '/' DK_normP(mid).name]);
                catch
                mid = find(ismember({DK_normP.name},[she_has.Var1{she_has.strainid==strainid},'.mat']));
                temp = load([DK_normP(mid).folder '/' DK_normP(mid).name]);
                end
                else
                temp = load([DK_normP(mid).folder '/' DK_normP(mid).name]);
            end
            if isfield(temp, 'medianNorm')
                normProfile= temp.medianNorm;
            elseif isfield(temp, 'normProfile')
                normProfile= temp.normProfile;
            end
        else
            temp = load([datafiles_rp(strains.bestrp(strainid)).folder '/' datafiles_rp(strains.bestrp(strainid)).name]);
            normProfile= temp.normProfile;
        end

        subplot(4,6,j-30)
        [out,sc4red]=mer_occupancy(normProfile(:,1),7,'window',21);
%         in out you have 8192 7mers and binding scores
        [~,top10]=maxk(out.score,10);

        [finalFreq]= seqLogoMin(out.mers.seq(top10),out.score(top10),'showFigs',true);
        title(temp.name)
        sgtitle(TF{i});
        
    end
    exportgraphics(gcf,[TF{i},'DDWTinvivomotifs.pdf'],'ContentType','vector')
end
           
%% fig.S3C OPN binding
%% strongly bound promoters
cursel = acol(denovoSid(1:6,1:9));
cursel = cursel(strains.max_corr(cursel)>0.9);
denovoRepro = cumsum(sum(strains.max_corr(denovoSid(:,1:9))>0.9,1));
cur_sp = StrainSumProm(:,cursel);
%% select current draw
denovoTFid = unique(strains.TFId(strains.type>30 & strains.type<40));

denovoSid = [];
for i = 1:6
    c=0;
    for j = denovoTFid'
        c=c+1;
        if j <27
            denovoSid(7,c) = find(strains.TFId==j & strains.type == 0);  
            denovoSid(8,c) = DBDid(j);
        else
            denovoSid(7,c) = find(strains.TFId==15 & strains.type == 0);  
            denovoSid(8,c) = DBDid(15);
        end
        denovoSid(i,c) = find(strains.TFId==j & strains.type == i+30);
    end
end

curDBD = StrainSumProm(:,denovoSid(8,1:9));
curTF = StrainSumProm(:,denovoSid(7,1:9));
log_sum = log2(cur_sp+700); % supress noise
zscore= (log_sum-mean(log_sum,1,'omitnan'))./std(log_sum,[],1,'omitnan');
z_log = zscore>5;
%% targets
opn = load('opnScore') ;
targetList = find(Protype<3&~isnan(opn.opnScore));
load('GP')
geneList = GP.gene_infoR64.name(targetList);
log_sPT = log_sum(targetList,:); 

%% opn scores & GC content
opnScore = opn.opnScore(targetList);
load('PromoterGCcontent.mat')
gctar = GCcontent(targetList);

%% cumsum binding ordered by opn
[~,orderTar] = sort(opnScore,'descend');
cumsums{1} = cur_sp(targetList(orderTar),:);
cumsums{2} = curDBD(targetList(orderTar),:);
cumsums{3} = curTF(targetList(orderTar),:);
x=[];
y=[];
e=[];
lo=[];
hi=[];
for i= 1:3
    cumsumSP = cumsum(cumsums{i});
    cumsumSP = cumsumSP./(cumsumSP(end,:));
    x{i} = (1:size(cumsumSP,1))';
    y{i} = (median(cumsumSP,2));
    e{i} =  std(cumsumSP,[],2);
    lo{i}  = movmean(y{i}-e{i},200); 
    hi{i} = movmean(y{i}+e{i},200);
    cumsumSP = [];
end
%% plot fig.S3C
figure

axes('Position',[0.2 0.3 0.6 0.4])
for i =[3,1,2]
    hp(i) = patch([x{i}; x{i}(end:-1:1); x{i}(1)], [lo{i}; hi{i}(end:-1:1); lo{i}(1)], 'r');
    hold on;
    hl(i) = line(x{i},movmean(y{i},200));
end
xticks([])
ylabel('Culmutative binding')

set(hp(1), 'facecolor', [134 147 230]/255, 'edgecolor', 'none','facealpha',0.7);
set(hp(2), 'facecolor', [.5 .5 .5], 'edgecolor', 'none','facealpha',0.7);
set(hp(3), 'facecolor', [147 176 152]/255, 'edgecolor', 'none','facealpha',0.7);
set(hl, 'color', [.5 .5 .5],'LineWidth',1);
set(gca,'XScale','log')
xlim([1 numel(orderTar)])
legend(hp,{'IDR-DBD','DBD','TF'})

axes('Position',[0.2 0.15 0.6 0.05])
xx = log10([1:1:numel(orderTar)]);
yy = ones(size(xx));
imagesc(xx,yy,opnScore(orderTar)')
set(gca,'colormap',brewermap(1000,'-RdBu'),'ytick',[],'xtick',[],'XScale','log')
colorbar
clim([-3 3])
xlabel('Promoters ordered by OPN')
exportgraphics(gcf,'DD_opn_log.eps','ContentType','vector')