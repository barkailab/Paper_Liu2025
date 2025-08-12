%% prepare 
% load('0724.mat','denovoinfo')
% get refsel and references
for i = 1:3
    refsel(i) = min(find(contains(denovoinfo.newname,['Ref',num2str(i)])));
    References{i} = denovoinfo.StrainName{refsel(i)};
end
%% Fig4A-left plot reference scaffold
ReferecesinShehas = find(ismember(FIrstBatch.newname,References(1:3)));
RefAAcom = AAcompmat(:,flipud(ReferecesinShehas));
xticknames = {'N','Q','S','T','G','DE','Hyd'};
figure

c=0;
for i = [aa2int('NQSTG'),21,aa2int(['F'])]
    c=c+1;
    axes('Position',[0.05+0.05*c,0.1,0.05,0.8])

    imagesc(ones(1,3),1:3,RefAAcom(i,:)',[0,max(RefAAcom(i,:))])
    set(gca, 'colormap',cmap{i},'XTick',1,'XTickLabel',xticknames{c},'ytick',[],'xticklabelrotation',0)
    nums =string(RefAAcom(i,:)');
    nums(strcmp(nums,"0"))=" ";
    yline([0.5:1:3-0.5],'w')

    if i==aa2int('F')
        nums = [FIrstBatch.hydaa(FIrstBatch.Order(1:3))];
       
    end
    text(repmat(0.6,1,3),1:3,nums)
    
end
%% Fig4A-right draw my design
allseq= fastaread('msn2denovo_seq.fasta');
FLcp = repmat([0.9608 0.6 0.6],20,1);
FLcp(aa2int('F'),:)=repmat([0 0 0],1,1);
FLcp(aa2int('L'),:)=repmat([1 1 1],1,1);
cursel = find(ismember(denovoinfo.newtype,'2Fdel')); % fig.S4B you select 1FtoL 
figure
c=0;
for i = 1:numel(cursel)
    c=c+1;
    subplot(numel(cursel),1,c)
    AAseq = aa2int(nt2aa(allseq(find(ismember({allseq.Header},denovoinfo.StrainName{cursel(i)}))).Sequence));


    imagesc(AAseq,[1,21])
    set(gca, 'ytick',[],'xtick',[])
    colormap(FLcp)
    % axis equal
    clim([1,21])
    yline(1.5,'k')
    ylim([0.5,1.5])
    xlim([0,460])
     title(strrep(denovoinfo.StrainName{cursel(i)},'_a',' '))
end
xticks(1:100:401)

%% select current draw
curIntAl = {'Fdel','Ldel'}%,'FF','LL'}%,'FtoL','Fto2L'}; %% figS2AB selects curIntAl 3~6
curIntAl = {'FtoL','Fto2L'};
drawsel2 = [];
for k = 1:numel(curIntAl)
    curInt = curIntAl{k};
    sel = find(contains(denovoinfo.newtype,curInt));
     sel = sel(strains.max_corr(denovoinfo.strainid(sel))>0.9);
    [~,idx] = sort(denovoinfo.denovoid(sel));

% get drawsel with the reference
    drawsel = [];
    for i = unique(denovoinfo.ref(sel))'
        if isempty(find(ismember(sel,refsel(i))))&k==1 
            cursel = [refsel(i);sel(denovoinfo.ref(sel)==i)];
            [a,b] = sort((cursel),'ascend');
            cursel = cursel(b);
            drawsel = [drawsel;cursel];
            cursel = [];
        else
            cursel = sel(denovoinfo.ref(sel)==i);
             [a,b] = sort((cursel),'ascend');
            cursel = cursel(b);
            drawsel = [drawsel;cursel];
            cursel = [];
        end
    end
    drawsel2=[drawsel2;drawsel];
end
drawsel = drawsel2;
%% Fig4B draw reproducibility
codefornFnL = arrayfun(@(x,y) ['F',num2str(x),'L',num2str(y)],denovoinfo.nF(drawsel),denovoinfo.nL(drawsel),'uniformoutput',false);
nn23 = [0    16    32    48    63    78    94   110   126];

tempmaxcorr = strains.max_corr(denovoinfo.strainid(drawsel));
tempmaxcorr(tempmaxcorr<0.9) = 0.5*tempmaxcorr(tempmaxcorr<0.9);
figure
scatter(codefornFnL, denovoinfo.ref(drawsel),  10+50*tempmaxcorr,'filled') 
xticks([sort(-nn23),nn23(2:end)])
xticklabels([sprintf('%d\n', sort(nn23,'descend')),sprintf('%d\n', nn23(2:end))])
ylim([1 3.5])
set(gca,'YDir','reverse')
hold on
scatter([1 1 1],[1 2 3],10+50*([1 0.4 0.1]))

%% Fig 4C Corr(most different F0) as Ref and corr F0
change = {{'Fdel','FF'},{'Ldel','LL'},'FtoL','Fto2L'};

clf
for i = 1:3
    for  j=1:2
    cursel = drawsel(denovoinfo.ref(drawsel)==i&contains(denovoinfo.newtype(drawsel),change{j}));

    cursel = cursel(strains.max_corr(denovoinfo.strainid(cursel))>0.9);
    if isempty(cursel)
        continue
    end
    id0 = find(contains(denovoinfo.newname,{'F0'})&denovoinfo.ref==i);
    refid = find(mod(denovoinfo.denovoid,1000)==0&denovoinfo.ref==i);
    refid = repmat(refid,numel(cursel),1);
    codefornFnL = denovoinfo.nF(cursel) + denovoinfo.nL(cursel);
    codefornFnL(contains(denovoinfo.newtype(cursel),'to')) = denovoinfo.nF(cursel(contains(denovoinfo.newtype(cursel),'to')));

    corrF0 = corr(StrainSumProm(Protype<3,unique(denovoinfo.strainid((id0)))),StrainSumProm(Protype<3,denovoinfo.strainid(cursel)),'rows','pairwise');
    
    seltemp= find(contains(denovoinfo.newtype(cursel),'Fdel'));
    [~,diffF01] = min(corrF0(contains(denovoinfo.newtype(cursel),'Fdel')));
    refid(contains(denovoinfo.newtype(cursel),'Fdel')) = cursel(seltemp(diffF01));
    
    seltemp= find(contains(denovoinfo.newtype(cursel),'Ldel'));
    [~,diffF02] = min(corrF0(contains(denovoinfo.newtype(cursel),'Ldel')));
    refid(contains(denovoinfo.newtype(cursel),'Ldel')) = cursel(seltemp(diffF02));

    corrDiffF0 = diag(corr(StrainSumProm(Protype<3,denovoinfo.strainid((refid))), ...
    StrainSumProm(Protype<3,denovoinfo.strainid(cursel)),'rows','pairwise'));

    if j==2
        colorid = '#FFDD99';
    elseif j ==1
        colorid = '#AED1B0';
    end

    subplot(1,3,i)
    
    hold on
    for kn = 1:numel(codefornFnL)
        scatter(codefornFnL(kn),corrF0(kn),125,denovoinfo.mutationid((cursel(kn))),'filled','MarkerFaceAlpha',codefornFnL(kn)./max(codefornFnL),'MarkerEdgeColor',[.75 .75 .75])
        text(codefornFnL(kn),corrF0(kn),denovoinfo.newname{(cursel(kn))})
    end
    clim([1 22])
    colormap(allcolors)
    
    end
end
%% Fig 4D 
% as above, just plot scatter(corrDiffF0, corrF0 ...)
%% Fig 4E corr all strains 
change = {'F','L'};

figure
for i= 1:3
        subplot(1,3,i)
        curself = drawsel(denovoinfo.ref(drawsel)==i&contains(denovoinfo.newname(drawsel),change{1}));
        tempbadf = find(strains.max_corr(denovoinfo.strainid(curself))<0.9);
         
        cursell= drawsel(denovoinfo.ref(drawsel)==i&contains(denovoinfo.newname(drawsel),change{2}));
        tempbadl =find(strains.max_corr(denovoinfo.strainid(drawsel))<0.9)+numel(curself); 
        corrmat = corr(StrainSumProm(Protype<3,denovoinfo.strainid([curself;cursell])));
        corrmat(tempbadf,:) = zeros(size(corrmat(tempbadf,:)));
        corrmat(:,tempbadf) = zeros(size(corrmat(:,tempbadf)));
        imagesc(corrmat,[0,1])
        set(gca,'ytick',1:numel([curself;cursell]),'yticklabel',denovoinfo.newname([curself;cursell]),...
            'xtick',1:numel([curself;cursell]),'xticklabel',denovoinfo.newname([curself;cursell]))
        % adjust colormap
        lowerlimwant = 0.7;
        colormap([.75 .75 .75;repmat([1 1 1],lowerlimwant*1000,1);brewermap(round((1-lowerlimwant)*1000),'Blues')])
        title(Refereces{i})
        axis equal
        ylim(xlim)
        xline([0.5:1:numel([curself;cursell])],'w')
        yline([0.5:1:numel([curself;cursell])],'w')
end