%% fig.S3A
TFdescription.TFi = cellfun(@(x) [x(1),lower(x(2:end))],TFdescription.TF,'UniformOutput',false);
colorblind_rgb = [
    0.0000 0.4500 0.7000   % dark blue
    0.8350 0.3690 0.0000   % vermillion
    0.0000 0.6190 0.4500   % bluish green
    0.9490 0.8940 0.1250   % yellow
    0.3500 0.3500 0.3500   % grey
    0.8000 0.4750 0.6550   % reddish purple
    0.3660 0.6000 0.8000   % sky blue
    0.9000 0.6000 0.0000   % orange
    0.9000 0.6000 0.6000   % pink
    0.8000 0.6000 0.4000   % light brown
];
figure
tempcmap = [166 207 227; 237 216 61;114 90 193;15 142 204;48 161 72;245 153 153]./255;
typei= 2;
    % subplot(1,2,typei)
    allval = strains.meanAbs(:,typei);
shiftScore = allval(drawlist(DDid(1:end-1),[31:36]));

 Goodscoremat = shiftScore./allval(drawlist(DDid(1:end-1),2));
 Goodscoremat(maxcorr(DDid(1:end-1),[31:36])<0.9)=-4;
[~,idx] = sort(max(Goodscoremat,[],2),'descend');
GoodxFLmat = corrWT(DDid(1:end-1),[31:36]);
GoodxFLmat(maxcorr(DDid(1:end-1),[31:36])<0.9)=0;
% corr. scaffold only
% Scaffoldid = [501,501,501,501,318,118];
Scaffoldid = drawlist(DDid(1:end-1),1);
GoodxScmat=[];
for i = 1:9
GoodxScmat(i,:) = corr(StrainSumProm(Protype<3, drawlist(DDid(i),31:36)), StrainSumProm(Protype<3,Scaffoldid(i)));
GoodxScmat(i,maxcorr(DDid(i),31:36)<0.9)=0;

end
GoodxFLmat = GoodxScmat;
hold on
n=0;
for i = idx'%[7,5,1,6,3,4,2,8,9]
    hold on
    n=n+1;
    scatter(n*ones(6,1),log(Goodscoremat((i),:)),100*GoodxFLmat((i),:)+1,1:6,'filled', ...
        'MarkerEdgeColor','#635A71','MarkerFaceAlpha',0.7)
    colormap(tempcmap)

    caxis([1 6])
    xlim([0 9])

end
xticks(1:9)
xticklabels(TFdescription.TFi(DDid(idx)))
ylabel('log(mutant/DBD-only promoter binding)')
title(names{typei})
  ylim([-1 1])
c=colorbar;
c.Ticks = 1+5/12:5/6:6;
c.TickLabels = allgrammarlist(31:36);
 scatter(9*ones(2,1),[0.8 1],[50,100],'k', ...
        'MarkerEdgeColor','#635A71')
 %% fig S3B
 %% Fig7 zscores (choose each TF)
tempcur_sp = StrainSumProm(:,drawlist(DDid(1:end-1),1));
templog_sum = log2(tempcur_sp+700); % supress noise
tempzscore= (templog_sum-mean(templog_sum,1,'omitnan'))./std(templog_sum,[],1,'omitnan');

n=0;
tfracstor=[];
tfMotifFracByz=[];
for tfi = 1:numel(DDid)-1
n=0;
for thr = 0:0.1:10
    n=n+1;
    % get TF targets
z_logDDFL = tempzscore>thr;
selTFsig = (z_logDDFL(:,tfi));
% how many denovo also bind them? check heree!!!
sel = find(strains.max_corr>0.9&ismember(strains.type,[31:36]));
cur_sp = StrainSumProm(:,sel);

log_sum = log2(cur_sp+700); % supress noise
zscore= (log_sum-mean(log_sum,1,'omitnan'))./std(log_sum,[],1,'omitnan');

z_DN_sig = zscore>thr;% DD > threshold
selDNsig = (sum(z_DN_sig,2)>1); % sel tar in DD
TF0DNsig = selTFsig.*selDNsig; % common targets
numTF0DNsig = sum(TF0DNsig(Protype<3,:))./sum(selDNsig(Protype<3)); 
DN0TFsig = selDNsig.*selTFsig;
numDN0TFsig = sum(DN0TFsig(Protype<3,:))./sum(selTFsig(Protype<3));

bothcaught= sum(selTFsig(Protype<3).*selDNsig(Protype<3));
tfracstor{tfi}(n,1) = thr;
tfracstor{tfi}(n,2) = bothcaught./sum(selTFsig(Protype<3));
tfracstor{tfi}(n,3) = sum(selTFsig(Protype<3));

% motif fraction
sigmotifYes = sum(motifYesMat(:,sel).*z_DN_sig(Protype<3,:))./sum(z_DN_sig(Protype<3,:));
tfMotifFracByz{tfi}(n,:) = sigmotifYes;

end
end
%% scatter of lines by each tf
figure

for i = 1:numel(DDid)-1
    subplot(2,5,i)
    hold on
    fill([0 0 3 3], [0 1 1 0],'g','EdgeColor','none','FaceColor','#cde6c6')

scatter(tfracstor{i}(tfracstor{i}(:,4)>0,1),tfracstor{i}(tfracstor{i}(:,4)>0,2),10,mean(tfMotifFracByz{i}(tfracstor{i}(:,4)>0,:),2),'filled','MarkerFaceAlpha',0.5)
hold on
colormap(brewermap(100,'Greens'))
clim([0 1])

xlabel('zscore')
ylabel('Fraction')%% old wrong
title(TFdescription.TF(DDid(i)))
xlim([0 5.5])
end
legend({'Targets threshold','Fraction of native targets'})
%% fig.S3C
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

%% fig. S3D seqlogo -FL
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
           
%% fig.S3E OPN binding
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
%% plot fig.S3E
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
%% Fig. S3G Euclidean distance from unbiased motifs and in-vitro motifs
tempsel = acol(drawlist(DDid(1:end-1),[31:36]));
DDselinDK = find(ismember([DK_normP.strainid], tempsel));

mybinedges = 0:0.02:1;
figure 
subplot(2,2,1)
histogram(UnbMCI(DDselinDK),'Normalization','probability','BinEdges',mybinedges,'FaceColor','#a8cedf','EdgeColor','none')

ylabel("Probability")
title('Denovo-DBD to in-vitro')
hold on
DselinDK = find(ismember([DK_normP.strainid], drawlist(DDid(1:end-1),2)));
subplot(2,2,2)
histogram([withref{DselinDK}],'normalization','probability','BinEdges',mybinedges,'FaceColor','#a8cedf','EdgeColor','none')

ylabel("Probability")
title('DBD-only to in-vitro')

subplot(2,2,3)
DFLselinDK = find(ismember([DK_normP.strainid], drawlist(DDid(1:end-1),1)));

histogram([withref{DFLselinDK}],'Normalization','probability','BinEdges',mybinedges,'FaceColor','#a8cedf','EdgeColor','none')

ylabel("Probability")
title('FL Euclidean to in-vitro')

% legend({'Denovo-DBD with invitro','DBD with invitro','FL with invitro'})
xlabel('Euclidean motif distance')
subplot(2,2,4)
histogram(cat(2,withother{DDselinDK}),'Normalization','probability','BinEdges',mybinedges,'FaceColor','#349946','FaceAlpha',0.6,'EdgeColor','none')

ylabel("Probability")
title('Denovo-DBD to other in vitro')
xlabel('Euclidean motif distance')

%% Fig. S3H: fraction of motif-containing targets
motifFrScore = MotifTarFrac(drawlist(DDid(1:end-1),[1:2,31:36]));
figure
hold on
% hold on
for i = 1:9
    hold on

    scatter(i*ones(8,1),1:8,20+100*adjmaxcorr(DDid(i),[1,2,31:36]),motifFrScore((i),:),'filled','MarkerEdgeColor','#635A71')
    colormap(brewermap(1000,'Greens'))

    clim([0.6 1])
    xlim([0 11])
    ylim([0 9])
    set(gca,'YDir','reverse')
    xticks([])
    yticks([])
    % pause()
end
axis off
