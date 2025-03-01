%% Get all orthologs for interested TFs
% %% download orthologs from YGOB
% % 
% TFinfo = readtable('TFdescription.xlsx');
% for i = 1:numel(TFinfo)
%     ygob{i} = getYGOB(TFinfo.TF{i});
%     ygob{i}.Header = strcat(ygob{i}.speciesPara,'_',ygob{i}.Header);
%     struc{i} = table2struct(ygob{i}(:,{'Header','Sequence'}));
%     fastawrite(fullfile(path,[TFinfo.TF{i},'_orthologs.fasta']),struc{i})
% end 
% here you can submit to clustral Omega for alignment
% https://www.ebi.ac.uk/jdispatcher/msa/clustalo

% save path
path = '/home/labs/barkailab/jingliu/Documents/Finalfiles/orthologs';
%%%% save all S.cer full length fasta
TF_AAseq = struct([]);
for i = 1:25
       seqfile = fastaread([path,TFinfo.TF{i},'_orthologs.fasta']);
    seqFL = seqfile(contains({seqfile.Header},'scer_A')).Sequence;
    TF_AAseq(i).Header = TFinfo.TF{i};
    TF_AAseq(i).Sequence = seqFL;
end
fastawrite('SelectedTFs_AAseq.fasta',TF_AAseq)
clearvars seqfile seqFL path
%% get lengths for Full length(FL), nonDBD, and DBD

load('MetaPredDisO.mat') % metapredict for disorder tendency
divyatf = readtable('descTF.xlsx');
for i = 1:25
    seqfile = fastaread([path,TFinfo.TF{i},'_orthologs.fasta']);
    seqFL = seqfile(contains({seqfile.Header},'scer_A')).Sequence;

    if TFinfo.nonDBD_start(i)<10
        nonDBDst = 1;
    else
        nonDBDst = TFinfo.nonDBD_start(i);
    end
    nonDBD = seqFL([nonDBDst:TFinfo.nonDBD_end(i)]);

    if TFinfo.FL_rm(i)>0
        DBDsel = true(1,length(seqFL));        
        DBDsel(sort([nonDBDst:TFinfo.nonDBD_end(i),TFinfo.FL_rm(i):TFinfo.FL_rmend(i)]))=false;
        DBD = seqFL(DBDsel);
        FLsel = true(1,length(seqFL));        
        FLsel(TFinfo.FL_rm(i):TFinfo.FL_rmend(i))=false;
        
        trueFL = seqFL(FLsel);
    else
        DBDsel = true(1,length(seqFL));        
        DBDsel(sort(nonDBDst:TFinfo.nonDBD_end(i)))=false;
        DBD = seqFL(DBDsel);
        tempid = find(ismember(divyatf.TF,TFinfo.TF{i}));
        if ~isempty(tempid)
            TFinfo.trDBD_st(i) = divyatf.FinalStart(tempid);
            TFinfo.trDBD_end(i) = divyatf.FinalEnd(tempid);
        end
        trueFL  = seqFL;
        FLsel = true(1,length(seqFL)); 
    end

    TFinfo.DBD_L(i) = length(DBD);
    TFinfo.trFL_L(i) = length(trueFL);
    TFinfo.nonDBD_L(i) = length(nonDBD);
    TFinfo.nonDBD{i} = nonDBD;
    TFinfo.disTMeta{i} = MetaPredDisO{i}([nonDBDst:TFinfo.nonDBD_end(i)]);
    TFinfo.disTMetaFL{i} = MetaPredDisO{i}(FLsel);
    TFinfo.nonDBD_OL(i) = TFinfo.nonDBD_L(i)-TFinfo.nonDBD_disL(i);    
end
clearvars seqfile seqFL nonDBDst DBDsel DBD Flsel trueFL 
%% Fig1B left show disorder tendency
[~,order] = sort(TFinfo.trFL_L(1:25),'descend');
UsedisO= TFinfo.disTMetaFL;
maxL = max(TFinfo.trFL_L);
figure
cn=0;
for i =order'
    cn=cn+1;
    subplot(26,1,cn)
    hold on
    c=[.75 .75 .75];
   if TFinfo.nonDBD_start(i)<10
       
        patch([TFinfo.nonDBD_end(i)+1 TFinfo.trFL_L(i)  TFinfo.trFL_L(i) TFinfo.nonDBD_end(i)+1],...
            [0 0 1 1] ,c,'EdgeColor','none')
    elseif TFinfo.nonDBD_start(i)>TFinfo.FL_rmend(i) & TFinfo.FL_rmend(i)>0
        rmd = TFinfo.FL_rmend(i)-TFinfo.FL_rm(i);
        patch([TFinfo.nonDBD_end(i)+1-rmd TFinfo.trFL_L(i)  TFinfo.trFL_L(i) TFinfo.nonDBD_end(i)+1-rmd],...
            [0 0 1 1],c,'EdgeColor','none')
     elseif TFinfo.nonDBD_end(i)<TFinfo.FL_rm(i)| TFinfo.nonDBD_start(i)>0
          patch([1 TFinfo.nonDBD_start(i)-1  TFinfo.nonDBD_start(i)-1 1],...
            [0 0 1 1],c,'EdgeColor','none')


    end
       patch([TFinfo.trDBD_st(i)  TFinfo.trDBD_end(i) TFinfo.trDBD_end(i) TFinfo.trDBD_st(i)],...
             [0 0 1 1],0.7*c,'EdgeColor','none')

      scatter(1:numel(UsedisO{i}),UsedisO{i},3,UsedisO{i}>0.5,'filled')
      yticks([])
      colormap([88,  89,  91; 247, 148, 29] / 255)
      xlim([0 maxL])
      ylabel(TFinfo.TF(i),'Rotation',0)
      ylim([0 1])
      xticks([])
end
xticks([0 800])
exportgraphics(gcf,'Fig1B_disT_orderL.pdf','ContentType','vector')
clearvars order UsedisO maxL cn i c rmd 
%% AA composition
aa2group = [8;1;5;2;8;5;2;7;1;4;4;1;8;3;7;6;6;3;3;4];

on_frq = [];
on_frqg=[];
for i = 1:25
    seq = aa2int(TFinfo.nonDBD{i});
     order_no = TFinfo.disTMeta{i}>0.5;
        order_yes = TFinfo.disTMeta{i}<0.5;
 
        on_frq(:,i) = accumarray(seq(order_no)',1,[20,1])/sum(order_no);
        on_frqg(:,i) = accumarray(aa2group,on_frq(:,i));
        oy_frq(:,i) = accumarray(seq(order_yes)',1,[20,1])/sum(order_yes);
        oy_frqg(:,i) =accumarray(aa2group,oy_frq(:,i));

        nonDBD_frq(:,i)= accumarray(seq',1,[20,1])/numel(seq);
        nonDBD_frqg(:,i) = accumarray(aa2group,nonDBD_frq(:,i));
end
clearvars seq order_no order_yes
%% get tf information in all yeast proteins
TFinfomation = readtable("TF_Information.txt");
UniqueTFDBid = unique(TFinfomation.DBID);
nonTFs = find(~ismember(GP.gene_infoR64.orf,UniqueTFDBid));
YesTFs = find(ismember(GP.gene_infoR64.orf,UniqueTFDBid));
%% generate nonTF disordered regions' AA fractions
load('seqOrdMeta.mat')
seq_noorder_noTF = {};
seq_noorder_yesTF = {};
Frq_noTF = zeros(20,1);
Frq_yesTF = zeros(20,1);
% -----for non-TFs-----
c=0;
for i = nonTFs'
    if numel(seqOrdMeta.Var1{i}(seqOrdMeta.Var2{i} > 0.5))>10
        c=c+1;
        seq_noorder_noTF{c} = seqOrdMeta.Var1{i}(seqOrdMeta.Var2{i} > 0.5);
        Frq_noTF(:,c) = Getfrq(aa2int(seqOrdMeta.Var1{i}(seqOrdMeta.Var2{i} > 0.5))',20);
        Frq_noTFg(:,c) = accumarray(aa2group,Frq_noTF(:,c));
        if Frq_noTF(:,c) == zeros(20,1)
            i %mark proteins without disordered regions
        end
    end
end

% -----get backgrounds-----
bg_nonTFs= aa2int(cat(2,seq_noorder_noTF{:}));
bg_yesTFs = aa2int(cat(2,seq_noorder_yesTF{:}));
bgg_noorder = aa2group(bg_noorder);
bgg_order = aa2group(bg_order);

Frq_nonTFs = Getfrq(bg_nonTFs',20);
% order_no_bgg = Getfrq(bgg_noorder,11);
Frq_yesTFs = Getfrq(bg_yesTFs',20);
order_yes_bgg = Getfrq(bgg_order,11);

allComp_g = accumarray(aa2group(1:20), allComp(1:20),[11,1]);
all_bg = allComp(1:20)/sum(allComp(1:20));
all_bgg = allComp_g/sum(allComp_g);

%% calculate AA enrichments
log2_on = log2(on_frq+0.001)-log2(median(Frq_noTF,2)+0.001);
log2_ong = log2(on_frqg+0.001)-log2(median(Frq_noTFg,2)+0.001);

clearvars allComp seqOrd seq_noorder seq_order bg_noorder bg_order bgg_noorder bgg_order allComp_g
%% AA groups color id
AAcolormap = repmat([180 192 192]/255,20,1);
AAcolormap(aa2int('KR'),:) = repmat([144 181 252]/255,2,1);
AAcolormap(aa2int('DE'),:) = repmat([227 98 151]/255,2,1);
AAcolormap(aa2int('FWY'),:) = repmat([209 250 208]/255,3,1);
AAcolormap(aa2int('LIV'),:) = repmat([247 255 135]/255,3,1);
AAcolormap(aa2int('STNQH'),:) = repmat([238 201 255]/255,5,1);
AAcolormap(21,:) = [1 1 1];
for i = 1:8
AAcolormapG(i,:) = AAcolormap(min(find(aa2group(1:20)==i)),:);
end

%% FigS1D overlook of AA fractions
figure
aaorder2 = aa2int(['FWY','LIV','NQ','ST','PG','MAC','KRH','DE']);
aa_names_g = {'KRH','DE','FWY','LIV','NQ','ST','PG','MAC'};

[~,order] = sort(TFinfo.trFL_L(1:25),'descend');
orderg = [3 4 5 6 7 8 1 2];
subplot(1,2,1)
imagesc(on_frq(aaorder2,order)',[0 0.2])
colormap(brewermap(1000,'Greens'))
xline(0.5:1:19.5,'w')
yline(0.5:1:24.5,'w')
set(gca,'xtick',1:20,'xticklabel',aa_names(aaorder2),'ytick',1:25,'yticklabel',TFinfo.TF(order),'XTickLabelRotation',0)
title('AA composition')
colorbar
subplot(1,2,2)
imagesc(on_frqg(orderg,order)',[0 0.3])
colormap(brewermap(1000,'Greens'))
set(gca,'xtick',1:9,'xticklabel',aa_names_g(orderg),'ytick',1:25,'yticklabel',TFinfo.TF(order))
title('AA composition')
xline(0.5:1:19.5,'w')
yline(0.5:1:24.5,'w')
colorbar
exportgraphics(gcf,'FigS1D_AAfrac.pdf','ContentType','vector')

%% Fig1B middle and FigS1D overlook of AA enrichments
clf
figure
subplot(1,2,1)
imagesc(log2_on(aaorder2,order)',[-4 4])
colormap(brewermap(1000,'-RdBu'))
xline(0.5:1:19.5,'w')
yline(0.5:1:24.5,'w')
set(gca,'xtick',1:20,'xticklabel',aa_names(aaorder2),'ytick',1:25,'yticklabel',TFinfo.TF(order))
title('AA composition')
colorbar
subplot(1,2,2)
imagesc(log2_ong(orderg,order)',[-2.5 2.5])
colormap(brewermap(1000,'-RdBu'))
set(gca,'xtick',1:9,'xticklabel',aa_names_g(orderg),'ytick',1:25,'yticklabel',TFinfo.TF(order))
title('AA composition')
xline(0.5:1:19.5,'w')
yline(0.5:1:24.5,'w')
colorbar
exportgraphics(gcf,'Fig1Bmid_AAenr.pdf','ContentType','vector')
%% Fig1B top show IDR background
figure
imagesc(median(Frq_noTFg(orderg,:),2)',[0 0.2])
colormap(brewermap(1000,'Greens'))
colorbar
title('All IDR protein composition(median)')

axis equal
ylim([0.5 1.5])
exportgraphics(gcf,'Fig1Btop_AAnonTFIDR.pdf','ContentType','vector')
%% select strains
drawlist_rp=[];
drawlist_ref =[];
drawlist = [];
clearvars abs_*
for i = 1:27
    a = find(ismember(strains.strain, TFdescription.FL{i}));
    b = DBDid(i);
    abs_log2{i,11} = 0;
    abs_dlog2{i,11} = 0;
    if ~isempty(strains.Abs_log2median(b))

        abs_log2{i,2} = strains.Abs_log2median(b);
        abs_dlog2{i,2} = strains.d_log2median(b);
        abs_short(i,2) = strains.d_log2median(b);
    end
    for j = [1:21,31:36]
        c = find(strains.type==j& strains.TFId ==i); %gramar repeats
        if isempty(c)
            continue
        end
        if numel(c)>1
            c = min(c);
            strains.strain(c)
        end
        drawlist_rp{i,j} = find(ismember(sTable.strainId ,c));
        if ~isempty(strains.Abs_log2median(c))
            abs_short(i,j) = strains.d_log2median(c);
        else
            abs_short(i,j) = 0;
        end
        if ~isempty(c)
            drawlist(i,j) = (c);
        end
    end
    drawlist_ref{i} = [a,b];
    drawlist(i,1:2) = [a,b];%% first col, WT, second, DBD,12 nonDBD, then as types
    nonDBDid = find(strains.type==12&strains.TFId==i);
    if ~isempty(nonDBDid)
        drawlist(i,12)=nonDBDid;
    end
end
drawlist(:,3:11)=0;
clearvars a b c 
%% get corrDBD WT for grammar
corrWT = zeros(27,36);
corrDBD =  zeros(27,36);
maxcorr =  zeros(27,36);
for i=[1:2,12:21,31:36]
    corrWT(drawlist(:,i)>0,i) = diag(corr(StrainSumProm(Protype<3,drawlist(drawlist(:,i)>0,1)),StrainSumProm(Protype<3,drawlist(drawlist(:,i)>0,i)),'rows','pairwise'));
    corrDBD(drawlist(:,i)>0,i) = diag(corr(StrainSumProm(Protype<3,drawlist(drawlist(:,i)>0,2)),StrainSumProm(Protype<3,drawlist(drawlist(:,i)>0,i)),'rows','pairwise'));
    maxcorr(drawlist(:,i)>0,i) = strains.max_corr(drawlist(drawlist(:,i)>0,i));
    strains.corrFL(drawlist(:,i)>0,i) = corrWT(drawlist(:,i)>0,i) ;
    strains.corrDBD(drawlist(:,i)>0,i) = corrDBD(drawlist(:,i)>0,i) ;
end
% for nonlinear dot sizes
adjmaxcorr = maxcorr;
adjmaxcorr(maxcorr<0.9) = 0.5*maxcorr(maxcorr<0.9);
%% allgrammarlist
allgrammarlist =[ 
    {'FL'            }
    {'DBD'           }
    {' '             }
    {' '             }
    {' '             }
    {' '             }
    {' '             }
    {' '             }
    {' '             }
    {' '             }
    {' '             }
    {' '             }
    {'LIVFWY toNQ'      }
    {'DE toNQ'      }
    {'LIV toF'       }
    {'FWY toL'       }
    {'KR toNQ'       }
    {'NQ toG'        }
    {'scramble'      }
    {'shift'         }
    {'scramble shift'}];
%% Fig1B right
figure
[~,order] = sort(TFinfo.trFL_L(1:25),'descend');
orderPert = [2 20 19 13 14 17 18 15 16];
allgrammarlist =strrep([{'FL','DBD'},repmat({' '},1,10),grammarlist],'_',' ');
hold on
for i = 1:9
    scatter(i*ones(25,1),[1:25]',adjmaxcorr(order,orderPert(i))*200+1,(corrWT(order,orderPert(i))),'filled','MarkerEdgeColor','k')
    clim([0 1])
end
colormap(brewermap(1000,'Purples'))
set(gca,'xtick',1:9,'xticklabel',allgrammarlist(orderPert),'ytick',1:25,'yticklabel',TFinfo.TF(order),'YDir','reverse','XTickLabelRotation',0)
colorbar
xlim([0 10])
ylim([1 28])
scatter([1,2,3],[26,26,26],[0.05 0.4 1]*200,[0 0 0],'filled')
text([1,2,3],[27,27,27],num2str([0.1 0.8 1]'))
%% output values in an excel
temptable = table();
temptable.TF = TFinfo.TF(1:25);
for i = 1:9
temptable.(allgrammarlist{orderPert(i)}) = corrWT(order,orderPert(i));
end
writetable(temptable,'/home/labs/barkailab/jingliu/Documents/Finalfiles/Sfiles/SuppValues.xlsx','Sheet','Figure1B')

%% correlation between scramble and shift
corrSS(find(drawlist(:,20)>0)) = diag(corr(StrainSumProm(Protype<3,drawlist(drawlist(:,20)>0,19)),StrainSumProm(Protype<3,drawlist(drawlist(:,20)>0,20)),'rows','pairwise'));

%% Fig1 DEF (FigS1 FG)
figure
n=0;
for i=[14,17,18,15,16,21,13]
    % get colorcode
    selAA = aa2int(extractBefore(allgrammarlist{i},' to'));
    temp_on = sum(on_frq(selAA,:));
    temp_bg = median(sum(Frq_noTF(selAA,:)),2);

    colorcode =log2(temp_on+0.001)-log2(temp_bg+0.001);
    n=n+1;
    subplot(2,4,n)
    hold on
    if i~=21
        scatter(corrDBD(maxcorr(:,i)>=0.9,i),corrWT(maxcorr(:,i)>=0.9,i),60,colorcode(maxcorr(:,i)>=0.9), ...
            "filled",'MarkerEdgeColor','k')
        scatter(corrDBD(20,i),corrWT(20,i),50*adjmaxcorr(20,i)+10,'markeredgecolor','r')
        text(corrDBD(maxcorr(:,i)>=0.9,i),corrWT(maxcorr(:,i)>=0.9,i),TF(maxcorr(:,i)>=0.9))
        set(gca,'Colormap',brewermap(1000,'-RdBu'),'CLim',[-2.5 2.5])
        if i == 13
            clim([-1 1])
            colorbar
        end
    % purple I used:4D2D84
    else
        scatter(corrDBD(maxcorr(:,i)>=0.9,i),corrWT(maxcorr(:,i)>=0.9,i),60, [0 0 0],"filled",'MarkerEdgeColor','none','MarkerFaceAlpha',0.5)%% no color
        text(corrDBD(maxcorr(:,i)>=0.9,i),corrWT(maxcorr(:,i)>=0.9,i),TF(maxcorr(:,i)>=0.9))
    end
 
    xlim([0,1])
    ylim([0,1])
    xticks([0 0.5 1])
    yticks([0 0.5 1])
%     
    plot(xlim,xlim,'Color',[.75,.75,.75])
    
    title(strrep(grammarlist{i-12},'_',' '))
    r=fill([0 0 1 0],[0 0 1 1], [205 230 196]/255);
    uistack(r,'bottom')
    xlabel('corr DBD')
    ylabel('corr FL')
end
sgtitle('Effect of DBD-only and pertubations (color: Swapped AA enrichment)')
% compare scramble and shift
n=n+1;
subplot(2,4,n)

scatter(corrWT(maxcorr(:,20)>=0.9,19),corrWT(maxcorr(:,20)>=0.9,20),60,corrSS(maxcorr(:,20)>=0.9),"filled",'MarkerEdgeColor','none')
hold on
scatter(corrWT(20,19),corrWT(20,20),50*maxcorr(20,21)+10,'markeredgecolor','r')
  
text(corrWT(maxcorr(:,20)>=0.9,19)+0.015,corrWT(maxcorr(:,20)>=0.9,20),TF(maxcorr(:,20)>=0.9))
  
title('Scrambel shift corr to respective FL (color: corr. scramble & shift)')
colormap(gca,brewermap(1000,'Purples'))
clim([0,1])
title('Compare scramble/shift effects')
xlabel('scramble',  'FontSize', 12);
 
xlim([0,1])
ylim([0,1])
xticks([0 0.5 1])
yticks([0 0.5 1])

plot(xlim,xlim,'Color',[.75 .75 .75])
ylabel('shift','fontsize',12);
colorbar


