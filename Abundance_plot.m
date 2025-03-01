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