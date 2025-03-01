sel = find(strains.max_corr(denovoinfo.strainid(1:116))>0.9);
cur_sp = StrainSumProm(:,denovoinfo.strainid(sel));

log_sum = log2(cur_sp+700); % supress noise
zscore= (log_sum-mean(log_sum,1,'omitnan'))./std(log_sum,[],1,'omitnan');
z_log = zscore>3;
sum(any(z_log,2)&Protype<3)
%%
[~,~,pertubationinfo.varid(1:21)]= unique(pertubationinfo.var(1:21));
%% Get correlation and polyfit R for each from to changes
% variables = {'nDE','nF','nF','nF','nH','nNQ','nP','nST','nG','nH','nP','nST','nH','nP'};

polyfitR = zeros(height(Protype),height(pertubationinfo));
corrR = zeros(height(Protype),height(pertubationinfo));
corrP = zeros(height(Protype),height(pertubationinfo));
fd = zeros(height(Protype),height(pertubationinfo));

typename=[];
for i = 1:height(pertubationinfo)-4
     cursel = find(ismember(denovoinfo.newtype(sel),pertubationinfo.fulltype{i}));
     pertubationinfo.cursels{i} = cursel;
   targettemp = 1:height(Protype);%find(any(z_log(:,cursel),2));
   for j = 1:height(Protype)
       try
        xn = denovoinfo.(pertubationinfo.var{i})(sel(cursel));
        usdedvar = pertubationinfo.var{i};
    catch
        xn = denovoinfo.nST(sel(cursel));
        usdedvar = 'nST';
       end
       xn = denovoinfo.(usdedvar)(pertubationinfo.ref(i))-xn;
            output = polyfit(log_sum(j,cursel)',xn,1);
            polyfitR(j,i) = output(1);
   end
   targets{i} = targettemp;
   [corrR(targettemp,i),corrP(targettemp,i)] = corr(log_sum(targettemp,cursel)',xn);
   fd(targettemp,i) = max(log_sum(targettemp,cursel),[],2)-min(log_sum(targettemp,cursel),[],2);
end

%% scatter of fold change and corr
   figure
   bgcolormap = brewermap(1001,'purples');
   for  i = 1:21
       subplot(5,5,i)
       scatter(fd(targets{i},i),corrR(targets{i},i),30,'b','filled','MarkerFaceAlpha',0.5)
       caxis([-2 2])
       colormap(gca,brewermap(1000,'-RdBu'))
       tempx = [0 7];
        tempy = [-1 1];
        if tempx(2)>2
        rectangle('Position',[1 0.75 tempx(2)-1 0.25],'EdgeColor','red','FaceColor','none')
        end
        rectangle('Position',[1 -1 tempx(2)-1 0.5],'EdgeColor','blue','FaceColor','none')
        bg=patch([tempx(1),tempx(2),tempx(2),tempx(1)], [-1 -1 1 1],bgcolormap(round(1000*(1-pertubationinfo.maxFL(i))),:),'FaceAlpha',1,'EdgeColor','none');
       uistack(bg,'bottom')
       xlim(tempx)
       ylim(tempy)
        title(pertubationinfo.fulltype{i})
   end

%% scrren out good regulated promoters
commontar = zeros(numel(Protype),21);
for i = 1:height(pertubationinfo)-6
    subtarU{i} = find(fd(:,i)>1 & corrR(:,i)>0.5 & any(z_log(:,pertubationinfo.cursels{i}),2));
    subtarD{i} = find(fd(:,i)>1 & corrR(:,i)<-0.75 & any(z_log(:,pertubationinfo.cursels{i}),2));
    commontar(subtarU{i},i)=corrR(subtarU{i},i);
    commontar(subtarD{i},i) =corrR(subtarD{i},i);
end
targetall = find(any(commontar'));


%% give DBD promoters not bound by the reference
logref(:,[1,3]) = log2(StrainSumProm(:,denovoinfo.strainid(refsel([1,3])))+700);
logMsn2nDBD(:,[1,2]) = log2(StrainSumProm(:,drawlist(15,[1:2]))+700);
logScaffold(:,[1,3]) = log2(StrainSumProm(:,denovoinfo.strainid(find(ismember(denovoinfo.newname,{'1F0','3F0'}))))+700);
zDBD =  (logMsn2nDBD-mean(logMsn2nDBD,1,'omitnan'))./std(logMsn2nDBD,[],1,'omitnan');
zRef =  (logref-mean(logref,1,'omitnan'))./std(logref,[],1,'omitnan');
zScaffold = (logScaffold-mean(logScaffold,1,'omitnan'))./std(logScaffold,[],1,'omitnan');
scaextar{1} = find(zScaffold(:,1)>3&zDBD(:,2)<2  &Protype<3);
scaextar{3} = find(zScaffold(:,3)>3&zDBD(:,2)<2  &Protype<3);

DBDextar = find(zDBD(:,2)>3 & zRef(:,1)<2&zRef(:,3)<2 &Protype<3);

reftar{1} = find(zRef(:,1)>3&zDBD(:,2)<2  &Protype<3);
reftar{3} = find(zRef(:,3)>3&zDBD(:,2)<2  &Protype<3);
id0([1,3]) = denovoinfo.strainid(find(ismember(denovoinfo.newname,{'1F0','3F0'})));

%% Fig 5B,C show binding to DBD exclusive targets by all 
clf

markerstouse = {'o',' ','^'};
toidColormap = [141, 11, 65;
    250,142,27;
    0,146,56;
    222, 76, 155;
    63,88,167;
    255,255,255]./255;
toidColormap = [brewermap(5,'Set2');1 1 1];
x=toidColormap(2,:);
toidColormap(2,:) =toidColormap(5,:);
toidColormap(5,:) = x;
refcolormap = [166 206 226;1 1 1;246 153 153]./255;
refcolorratio = 1;
subplot(1,2,1)
hold on
xline(0,'--')
yline(0,'--')
for i = 1:21
    tempref = pertubationinfo.ref(i);
    bindDBD = median(log_sum(DBDextar,pertubationinfo.cursels{i}(end)));
    xn = pertubationinfo.xn{i};
     DBDfd = bindDBD-median(logref(DBDextar,tempref));
     xDBDenr = pertubationinfo.end_DBD(i);
%   Allfd = median(fd(subtarU,i));
     scatter(DBDfd,xDBDenr,100,pertubationinfo.varid(i),'filled','MarkerEdgeColor', ...
         toidColormap(pertubationinfo.toid(i),:),'LineWidth',1,'Marker',markerstouse(pertubationinfo.ref(i)))
    colormap(toidColormap)
    caxis([1 6])
%     text(DBDfd,xDBDenr,[pertubationinfo.type{i}])
 
end
hold on
for i = [1,3]
        bindTarDBD_ref = median(logref(DBDextar,i));
        bindTarDBD_id0 = median(logScaffold(DBDextar,i));
        id0xDBD = corr(StrainSumProm(Protype<3,drawlist(15,2)),StrainSumProm(Protype<3,id0(i)));
           scatter(bindTarDBD_id0-median(logref(DBDextar,i)),(id0xDBD),120,'k',markerstouse{i})
        scatter(bindTarDBD_ref-median(logref(DBDextar,i)),refDBD(i),120,'k',markerstouse{i},'filled')
end
xlabel('DBD targets enrichment')
ylabel('corrDBD')

subplot(1,2,2)
hold on

for i = 1:21
    tempref = pertubationinfo.ref(i);
    bindref = median(log_sum(reftar{tempref},pertubationinfo.cursels{i}(end)));
    xn = pertubationinfo.xn{i};
     reffd = bindref-median(logref(reftar{tempref},tempref));
     xrefenr = pertubationinfo.end_ref(i);
%   Allfd = median(fd(subtarU,i));
     scatter(reffd,xrefenr,100,pertubationinfo.varid(i),'filled','MarkerEdgeColor', ...
         toidColormap(pertubationinfo.toid(i),:),'LineWidth',1,'Marker',markerstouse(pertubationinfo.ref(i)))
    colormap(toidColormap)
    caxis([1 6])
%     text(reffd,xrefenr,[pertubationinfo.type{i}])
 
end
hold on
for i = [1,3]
        bindTarref_ref = median(logref(reftar{i},i));
        bindTarref_id0 = median(logScaffold(reftar{i},i));
        id0xref = corr(StrainSumProm(Protype<3,denovoinfo.strainid(refsel(i))),StrainSumProm(Protype<3,id0(i)));
           scatter(bindTarref_id0-median(logref(reftar{i},i)),id0xref,120,'k',markerstouse{i})
        scatter(0,1,120,'k',markerstouse{i},'filled')
end
xline(0,'--')
yline(0,'--')
xlabel('Ref targets enrichment')
ylabel('corrRef')
c=colorbar;
c.Ticks=1.5:0.8:5.5;
c.TickLabels = {'NQ','ST','G','H','P','DE'};
%% Fig 5E scatter of fold change and corr
figure

for  i = 1:21
    subplot(5,5,i)
    scatter(fd(targets{i},i),corrR(targets{i},i),30,'b','filled','MarkerFaceAlpha',0.5)
    clim([-2 2])
    tempx = [0 7];
    tempy = [-1 1];
    if tempx(2)>2
        rectangle('Position',[1 0.75 tempx(2)-1 0.25],'EdgeColor','red','FaceColor','none')
    end
    rectangle('Position',[1 -1 tempx(2)-1 0.5],'EdgeColor','blue','FaceColor','none')
    xlim(tempx)
    ylim(tempy)
    title(pertubationinfo.fulltype{i})
end
exportgraphics(gcf,'opnfdandr.pdf','ContentType','vector')
%% Fig 5F show examples in the main figure

for i =1:21
    clf
    cursel = find(ismember(denovoinfo.newtype(sel),pertubationinfo.fulltype{i}));
    cursel = [refsel(denovoinfo.ref(sel(cursel(1))));sel(cursel)];
    
    cur_sp = StrainSumProm(:,denovoinfo.strainid(cursel));
    log_sumT = log2(cur_sp+700); % log for msn2 & my pertubations
    [~,orderU] = sort(log_sumT(subtarU{i},end),'descend');
    [~,orderD] = sort(log_sumT(subtarD{i},1),'descend');
   
    try
        nHyd = denovoinfo.(pertubationinfo.var{i})((cursel));
        usdedvar = pertubationinfo.var{i};
    catch
        nHyd = denovoinfo.nST((cursel));
        usdedvar = 'nST';
    end
    nHyd = nHyd(1) - nHyd;
    clf
    axes('position',[0.05 0.47 0.14 0.4])
    targetList = subtarU{i}(orderU);
    imagesc(log_sumT(targetList,:),[12,20])
    
    hold on
    geneids1 = GP.gene_infoR64.name(targetList);
    set(gca,'XTick',[],'XTickLabel', denovoinfo.newname((cursel)),...
    'Ytick',1:numel(geneids1),'YTicklabel',geneids1,'XTickLabelRotation',90,'YDir', ...
    'reverse','ticklabelinterpreter','none','colormap',brewermap(1000,'Greens'),'fontsize',4)
    ylabel(sprintf('n = %d',numel(targetList)),'Rotation',0)
    title(pertubationinfo.fulltype{i})
    colorbar

    axes('position',[0.05 0.05 0.14 0.4])
    targetList = subtarD{i}(orderD);
    imagesc(log_sumT(targetList,:),[12,20])
    hold on
    geneids1 = GP.gene_infoR64.name(targetList);
    set(gca,'XTick',1:numel(cursel),'XTickLabel', denovoinfo.step((cursel)),...
    'Ytick',1:numel(geneids1),'YTicklabel',geneids1,'XTickLabelRotation',0,'YDir', ...
    'reverse','ticklabelinterpreter','none','colormap',brewermap(1000,'Greens'),'fontsize',4)
    colorbar

    ylabel(sprintf('n = %d',numel(targetList)),'Rotation',0)

    axes('Position',[0.3,0.15 0.17 0.7])

    hold on
    mean_bind = mean(log_sumT(subtarU{i},:),1);
    std_bind = std(log_sumT(subtarD{i},:));
    scatter(nHyd,mean_bind,10+nHyd,'b','filled')
    % errorbar(nHyd,mean_bind,std_bind)
    plot(nHyd(1:numel(cursel)),mean_bind(1:numel(cursel)),'Color','b')
    % yline(mean(logDBD(subtarU{i},5)),'r')
    
    title('top Pos')
    %  set(gca,'xdir','reverse')
    axes('Position',[0.5,0.15 0.17 0.7])
    
    hold on
    mean_bind = mean(log_sumT(subtarD{i},:));
    std_bind = std(log_sumT(subtarD{i},:));
    scatter(nHyd,mean_bind,10+nHyd,'b','filled')
    % errorbar(nHyd,mean_bind,std_bind)
    plot(nHyd(1:numel(cursel)),mean_bind(1:numel(cursel)),'Color','b')
    
    ylabel('log2 binding')
    title(['tail Pos'])

    exportgraphics(gcf,[pertubationinfo.fulltype{i},'promoters.pdf'],'ContentType','vector')

end
%% gene list
RegList = table();
for i = 1:21
    RegList = table();
    
    cursel = find(ismember(denovoinfo.newtype(sel),pertubationinfo.fulltype{i}));
    cursel = [refsel(denovoinfo.ref(sel(cursel(1))));sel(cursel)];
    
    cur_sp = StrainSumProm(:,denovoinfo.strainid(cursel));
    log_sumT = log2(cur_sp+700); % log for msn2 & my pertubations
    [~,orderU] = sort(log_sumT(subtarU{i},end),'descend');
    [~,orderD] = sort(log_sumT(subtarD{i},1),'descend');

    try
        nHyd = denovoinfo.(pertubationinfo.var{i})((cursel));
        usdedvar = pertubationinfo.var{i};
    catch
        nHyd = denovoinfo.nST((cursel));
        usdedvar = 'nST';
    end
    nHyd = nHyd(1) - nHyd;
    clf
    
    targetList = subtarU{i}(orderU);
    RegList.geneNames = GP.gene_infoR64.name(targetList);
    RegList.foldchange = fd(targetList,i);
    RegList.graduality = corrR(targetList,i);
    for j = 1:numel(cursel)
    RegList.(['log2Binding',num2str(j)]) = log_sumT(targetList,j);
    end
    writetable(RegList,'ScaffoldUpgenes.xlsx','Sheet',pertubationinfo.fulltype{i})
end
%% colormap for Fig6G
% Define the RGB values for baby blue, white, and pastel pink
baby_blue = [173, 216, 230] / 255;  % #ADD8E6
white = [1, 1, 1];                  % White
pastel_pink = [255, 192, 203] / 255; % #FFC0CB

% Number of colors in the colormap
num_colors = 256;

% Define positions of key colors
x = [1, num_colors/2, num_colors]; % Positions of defined colors
colors = [baby_blue; white; pastel_pink]; % Color matrix

% Interpolate for smooth gradient
r = interp1(x, colors(:,1), linspace(1, num_colors, num_colors));
g = interp1(x, colors(:,2), linspace(1, num_colors, num_colors));
b = interp1(x, colors(:,3), linspace(1, num_colors, num_colors));

% Combine into colormap
pastel_colormap = [r', g', b'];

% Apply colormap

%% Fig 6G , S5D plot the ones with much up /down
selperb = find(cellfun(@numel, subtarD)>10);
upordown = subtarD; %choose subtarD for down
[~,orderPerb] = sort(cellfun(@numel, upordown(selperb)),'descend');
orderPerb = selperb(orderPerb);
% order each U genes 
for i = 1:21
    [~,order] = sort(corrR(upordown{i},i),'descend');
upordown{i} = upordown{i}(order);
end

figure
c=0;
targetUtar = [];
commontarUPart = [];
plotPerb = orderPerb(find(cellfun(@numel, upordown(orderPerb))>10));
xv = 0.05;
totalgenetoshow = sum(cellfun(@(x) numel(x), upordown));
gap = 0.005;
genewidth = (1-gap*numel(order))./totalgenetoshow;

for k=1:numel(plotPerb)
    i=plotPerb(k);
    c=c+1;
    if k>1
    xv= xv+numel(upordown{plotPerb(k-1)})*genewidth+gap;
    end
    yv = 0.95;
    for j = orderPerb
      
       yv = yv-0.04-gap;
       axes('position',[xv yv numel(upordown{i})*genewidth 0.04])
       imagesc(commontar(upordown{i},j)',[-1 1])
       if i==j
           colormap(gca,brewermap(1000,'-RdBu'))
       else
           colormap(gca,pastel_colormap)
       end

       if i==orderPerb(1)
           ylabel(pertubationinfo.fulltype{j},'Rotation',0)
       end

       set(gca,'xtick',[],'ytick',[])
    end
    
end
exportgraphics(gcf,'/home/labs/barkailab/jingliu/Documents/Finalfiles/Sfiles/AllDownOverlap3.pdf','ContentType','vector')

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
