%% Fig S4C absolute measurements
change = {{'Fdel','FF'},{'Ldel','LL'},'FtoL','Fto2L'};
load('0724','allcolors')
clf
for typei =2
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

    corrF0 = strains.meanAbs(denovoinfo.strainid(cursel),typei);
    corrF0 = log(corrF0);
    strains.strain(denovoinfo.strainid(cursel))
   
    if j==2
        colorid = '#FFDD99';
           
    elseif j ==1
        colorid = '#AED1B0';
        
    end

    subplot(2,3,i+3*(j-1))
    
    hold on
    for kn = 1:numel(codefornFnL)
        scatter(codefornFnL(kn),corrF0(kn),125,denovoinfo.mutationid((cursel(kn))),'filled','MarkerFaceAlpha',codefornFnL(kn)./max(codefornFnL),'MarkerEdgeColor',[.75 .75 .75])
        % text(codefornFnL(kn),corrF0(kn),denovoinfo.newname{(cursel(kn))})
    end
    clim([1 22])
    colormap(allcolors)
  ylim([15.4 16.04])
  xlim([0 150])
    end
end
end
%% Fig S4D
figure
tempmaxcorr = MotifTarFrac(denovoinfo.strainid(drawsel));
codefornFnL = arrayfun(@(x,y) x-y,denovoinfo.nF(drawsel),denovoinfo.nL(drawsel),'uniformoutput',true);
xpos = 1:9;
scatter(codefornFnL, denovoinfo.ref(drawsel),150,tempmaxcorr,'filled') 
colormap(brewermap(1000,'Greens'))

clim([0.3 1])

xticks([sort(-xpos),xpos(2:end)])
xticklabels([sprintf('%d\n', sort(xpos,'descend')),sprintf('%d\n', xpos(2:end))])
ylim([1 3.5])
set(gca,'YDir','reverse')
