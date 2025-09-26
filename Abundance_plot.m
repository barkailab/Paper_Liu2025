%% abundance plot fig.S1E, S3D
selpert =[31:36];
absshow= abs_short(:,selpert);

colormapp = brewermap(1001,'Greens');
absshowcolor = round(1000*(absshow-min(absshow(:)))./(max(absshow(:))-min(absshow(:))));
figure
hold on
cn=0;
for i = find(any(absshow,2))'
    cn=cn+1;
    for j =1:numel(selpert)
        if absshow(i,j)~=0
            scatter(cn,j,20+100*adjmaxcorr(i,selpert(j)),absshowcolor(i,j)+1,'filled','MarkerEdgeColor','k')
        end
    end
end

set(gca,'ytick',1:numel(selpert),'YTickLabel',pertubationlist(selpert), ...
    'xtick',1:cn,'xticklabel',TF(find(any(absshow,2))),'YDir','reverse','colormap',colormapp)
c=colorbar;
xlim([-0.5,cn+1])
ylim([-0.5,numel(selpert)+1])
axis equal

%% abundance for S
figure
subplot(2,3,1)

selstr = find(ismember(strains.type,[13:20])&strains.Abs_log2median~=0&strains.TFId~=1);
xv= strains.max_corr(selstr);
yv = strains.d_log2median(selstr);

hold on
scatter(xv,yv,40,[0.2588    0.5725    0.7765],'filled')
[c,p] = corr(xv,yv);
[x,y] = polyfit(xv,yv,1);
x_fit = linspace(0,1,100);
y_fit = polyval(x, x_fit);

colormap(brewermap(1000,'Blues'))
caxis([0 1])

xlabel('Corr. Reference')
ylabel('Abandance')
title(sprintf('S1E: c = %.2f, p = %.2f >0.05*', round(c,1), round(p,1)))
subplot(2,3,2)

selstr = find(ismember(strains.type,[31:36])&strains.Abs_log2median~=0&strains.TFId~=1);
xv= strains.max_corr(selstr);
yv = strains.d_log2median(selstr);

hold on
scatter(xv,yv,40,[0.2588    0.5725    0.7765],'filled')
[c,p] = corr(xv,yv);
[x,y] = polyfit(xv,yv,1);
x_fit = linspace(0,1,100);
y_fit = polyval(x, x_fit);


colormap(brewermap(1000,'Blues'))
caxis([0 1])

xlabel('Corr. Reference')
ylabel('Abandance')
title(sprintf('S3F: c = %.2f, p = %.2f >0.05*', round(c,1), round(p,1)))

% abundance for S std
% figure
subplot(2,3,4)

selstr = find(ismember(strains.type,[13:20])&strains.Abs_log2median~=0&strains.TFId~=1);
xv= strains.max_corr(selstr);
yv = cellfun(@mean,strains.Abs_std(selstr));

hold on
scatter(xv,yv,40,[0.2588    0.5725    0.7765],'filled')
[c,p] = corr(xv,yv);
[x,y] = polyfit(xv,yv,1);
x_fit = linspace(0,1,100);
y_fit = polyval(x, x_fit);

colormap(brewermap(1000,'Blues'))
caxis([0 1])

xlabel('Corr. Reference')
ylabel('Abandance')
title(sprintf('S1E: c = %.2f, p = %.2f >0.05*', round(c,1), round(p,1)))
subplot(2,3,5)

selstr = find(ismember(strains.type,[31:36])&strains.Abs_log2median~=0&strains.TFId~=1);
xv= strains.max_corr(selstr);
yv = cellfun(@mean,strains.Abs_std(selstr));

hold on
scatter(xv,yv,40,[0.2588    0.5725    0.7765],'filled')
[c,p] = corr(xv,yv);
[x,y] = polyfit(xv,yv,1);
x_fit = linspace(0,1,100);
y_fit = polyval(x, x_fit);


colormap(brewermap(1000,'Blues'))
caxis([0 1])

xlabel('Corr. Reference')
ylabel('Abandance')
title(sprintf('S3F: c = %.2f, p = %.2f >0.05*', round(c,1), round(p,1)))

subplot(2,3,6)

selstr = find(ismember(strains.type,30)&strains.Abs_log2median~=0&strains.TFId~=1);
xv= strains.max_corr(selstr);
yv = cellfun(@mean,strains.Abs_std(selstr));

hold on
scatter(xv,yv,40,[0.2588    0.5725    0.7765],'filled')
[c,p] = corr(xv,yv);
[x,y] = polyfit(xv,yv,1);
x_fit = linspace(0,1,100);
y_fit = polyval(x, x_fit);


colormap(brewermap(1000,'Blues'))
caxis([0 1])

xlabel('Corr. Reference')
ylabel('Abandance')
title(sprintf('S5C: c = %.2f, p = %.2f >0.05*', round(c,1), round(p,1)))


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


