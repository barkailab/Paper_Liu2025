%% get meansurements
c=0;
for index = 1:height(datafiles_rp)
    c=c+1;
 
temp = load([datafiles_rp(index).folder '/' datafiles_rp(index).name]);
if ~isfield(temp,'name')
    temp.name = strrep(datafiles_rp(index).name,'.mat','');
end
    [~, AbsMat, TF] = Abs_binding(temp);
    if ~isempty(AbsMat)
        datafiles_rp(index).AbsMotif = AbsMat(1);
        datafiles_rp(index).AbsPromoter = AbsMat(2);
        datafiles_rp(index).AbsMitochondria = AbsMat(3);
        datafiles_rp(index).AbsTelomeric = AbsMat(4);
        datafiles_rp(index).TFuse = TF;
    end

end
%% correct shape 
for i = 1:height(datafiles_rp)
    if size(datafiles_rp(i).MeanAbs,2)~=1
        datafiles_rp(i).MeanAbs = datafiles_rp(i).MeanAbs';
    end
end
%% show binding by each strain
for i = 1:height(strains)
    % if ~ismember(strains.people{i},'divyakr')
        strains.AbsSel(i)=1;
        tempairs = strains.pairs{i};
        temprepeats = unique(tempairs(:));
    if numel(tempairs)==2 & strains.max_corr(i)>=0.9
        withGoodReads = sTable.reads(tempairs)>1e6;
        strains.GoodReadsCorr{i} = strains.all_corr{i}(withGoodReads);
       goodrepeats = temprepeats(sTable.reads(temprepeats)>=1e6);
        badrepeats = temprepeats(sTable.reads(temprepeats)<1e6);
        strains.goodAbsSum{i} = [[datafiles_rp(goodrepeats).AbsMotif];[datafiles_rp(goodrepeats).AbsPromoter];[datafiles_rp(goodrepeats).AbsMitochondria];[datafiles_rp(goodrepeats).AbsTelomeric]]';
    
        strains.badAbsSum{i} = [[datafiles_rp(badrepeats).AbsMotif];[datafiles_rp(badrepeats).AbsPromoter];[datafiles_rp(badrepeats).AbsMitochondria];[datafiles_rp(badrepeats).AbsTelomeric]]';
        temp = [datafiles_rp(goodrepeats).MeanAbs]';
        if ~isempty(temp)
            strains.goodAbsMean{i} = temp;
        end
        temp = [datafiles_rp(badrepeats).MeanAbs]';
        if ~isempty(temp)
            strains.badAbsMean{i} = temp;
        end
        
    else
    withGoodReads = sum(sTable.reads(tempairs)>1e6,2)==2;
    if ~isempty(withGoodReads)
    strains.GoodReadsCorr{i} = strains.all_corr{i}(1,withGoodReads);
    end
    % goodrepeats = unique(tempairs(withGoodReads,:));
    % badrepeats = unique(tempairs(~ismember(tempairs,goodrepeats)));
      goodrepeats = temprepeats(sTable.reads(temprepeats)>=1e6);
        badrepeats = temprepeats(sTable.reads(temprepeats)<1e6);
      
    strains.goodAbsSum{i} = [[datafiles_rp(goodrepeats).AbsMotif];[datafiles_rp(goodrepeats).AbsPromoter];[datafiles_rp(goodrepeats).AbsMitochondria];[datafiles_rp(goodrepeats).AbsTelomeric]]';

    strains.badAbsSum{i} = [[datafiles_rp(badrepeats).AbsMotif];[datafiles_rp(badrepeats).AbsPromoter];[datafiles_rp(badrepeats).AbsMitochondria];[datafiles_rp(badrepeats).AbsTelomeric]]';
        temp = [datafiles_rp(goodrepeats).MeanAbs]';
        if ~isempty(temp)
            strains.goodAbsMean{i} = temp;
        end
        temp = [datafiles_rp(badrepeats).MeanAbs]';
        if ~isempty(temp)
            strains.badAbsMean{i} = temp;
        end
        temp=[];
    end
       % end
end
yline(0.9)
%% look at all mean & std of abs binding
stdAll = zeros(705,4);
stdGood = zeros(705,4);
stdBad = zeros(705,4);
MAll = zeros(705,4);
MGood = zeros(705,4);
MBad = zeros(705,4);
for i = 1:height(strains)
    if strains.AbsSel(i)==1
        [stdAll(i,:),MAll(i,:)] = std([strains.goodAbsSum{i};strains.badAbsSum{i}]);
        if height(strains.goodAbsSum{i}>2)
        [stdGood(i,:),MGood(i,:) ] = std(strains.goodAbsSum{i});
       
        end
        if height(strains.badAbsSum{i}>2)
        [stdBad(i,:),MBad(i,:)] =std(strains.badAbsSum{i});

        end
    end
end
strains.meanAbs = MAll;

%% full gradient changes
for i = 1:height(pertubationinfo)
    PerbMat(i,1:numel(pertubationinfo.cursels{i})) = strains.meanAbs(denovoinfo.strainid(sel(pertubationinfo.cursels{i})),2);   
end