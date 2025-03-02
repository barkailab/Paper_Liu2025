function organiseDataGUI()
    h.fig=figure;
    h.repeatFolder=uicontrol('Style','pushbutton','String','select Folder for sample profiles','Units','normalized','Position',[0.05 0.9 0.25 0.03],'Callback',@selRepeatFolder);
    h.importOutfolder=uicontrol('Style','pushbutton','String','Import more outFiles','Units','normalized','Position',[0.05 0.85 0.25 0.03],'Callback',@importOut);
    h.expList=uicontrol('Style','listbox','String',{},'Units','normalized','Position',[0.35 0.1 0.35 .825]);
    h.meanFolder=uicontrol('Style','pushbutton','String','select Folder for mean profiles','Units','normalized','Position',[0.73 0.9 0.25 0.03],'Callback',@selMeanFolder);
    h.createMeans=uicontrol('Style','text','String','available experiments (# Repeats / mean status)','Units','normalized','Position',[0.35 0.93 0.35 0.025]);
    h.checkRpts=uicontrol('Style','pushbutton','String','Check Repeats','Units','normalized','Position',[0.73 0.85 0.25 0.03],'Callback',@checkRepeats);    
    h.addMeans=uicontrol('Style','pushbutton','String','Add mean profiles','Units','normalized','Position',[0.73 0.80 0.25 0.03],'Callback',@(a,b)calcMeans(a,b,'add'));
    h.allMeans=uicontrol('Style','pushbutton','String','recalculate all Means','Units','normalized','Position',[0.73 0.75 0.25 0.03],'Callback',@(a,b)calcMeans(a,b,'all'));
    h.allMeans=uicontrol('Style','pushbutton','String','calcualte selected Means','Units','normalized','Position',[0.73 0.70 0.25 0.03],'Callback',@(a,b)calcMeans(a,b,'selection'));
    
    guidata(h.fig,h)
end

function selRepeatFolder(a,b)
    repeatFolder=uigetdir('~');
    h=guidata(a.Parent);
    h.repeatFolder.String=repeatFolder;
    guidata(h.fig,h)
    refreshExpList(h,[])
end

function importOut(a,b)
    outFolder=uigetdir('/home/labs/barkailab/LAB/data/SEQ/')
    if exist([outFolder,'/WellList.xlsx'],'file')
        if numel(dir([outFolder,'/*stat*.tsv']))==1
            [fullProfile,sample,smeta]=import_outfolder(outFolder,'sample','WellList','meta',true,'sort',false);
        else
            [fullProfile,sample,smeta]=import_outfolder(outFolder,'sample','WellList','meta',false,'sort',false);
        end
    elseif exist([outFolder,'/AdaptaramaSampleSheet.xlsx'],'file')
         if numel(dir([outFolder,'/*stat*.tsv']))==1
            [fullProfile,sample,smeta]=import_outfolder(outFolder,'sample','AdaptaramaSampleSheet','meta',true,'sort',false);
        else
            [fullProfile,sample,smeta]=import_outfolder(outFolder,'sample','WellList','meta',false,'sort',false);
        end       
    else
        indexFiles=dir([outFolder,'/*.xlsx'])
        if numel(indexFiles)>0
             [fullProfile,sample,smeta]=import_outfolder(outFolder,'sample','Gilad','meta',true,'sort',false);
        else
            [fullProfile,sample,smeta]=import_outfolder(outFolder,'meta',true,'sort',false);
        end
    end
    GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
    gb=true(size(fullProfile,1),1);
    gb(GP.chrIdx(17)+1:GP.chrIdx(18))=false;
    for i=1:size(GP.varRegions)
        posI=sort([GP.chrIdx(GP.varRegions{i,1})+GP.varRegions{i,2:3}]);
        gb(posI(1):posI(2),:)=false;
    end
    totalReads=sum(fullProfile(gb,:));
    normProfile=fullProfile./totalReads*sum(gb);
    clear fullProfile gb posI i 
    load('/home/labs/barkailab/tamarge/Master/mat Files/promoterLengthsORF.mat','promoterLengthsORF');
    metaPro=metaProfilePromLenDivya(normProfile,'promEnd','position','afterTss',0,'promLen',promoterLengthsORF);
    sumProm=squeeze(sum(metaPro,2,'omitnan'));
    sumProm(isnan(promoterLengthsORF),:)=NaN;
    clear temp
    h=guidata(a.Parent)
    for i=1:numel(sample)
        temp.normProfile=normProfile(:,i);
        temp.sumProm=sumProm(:,i);
        temp.totalReads=totalReads(i);
        temp.name=sample{i};
        temp.meta=smeta(i,:);
        temp.file=[outFolder ,'/', temp.meta.File{1}]
        save(sprintf('%s/%s.mat',h.repeatFolder.String ,temp.name),'-struct','temp','-v7.3')
        clear temp
    end
    guidata(h.fig,h);
    refreshExpList(h,unique(regexp(sample,'.*(?=_)','match','once')))
end

function selMeanFolder(a,b)
    meanFolder=uigetdir('~');
    h=guidata(a.Parent);
    h.meanFolder.String=meanFolder;
    guidata(h.fig,h)
    refreshExpList(h,[])
end

function refreshExpList(h,selExp)
    allRpts=struct2table(dir([h.repeatFolder.String,'/*.mat']));
    if size(allRpts,1)>0
        [exps,~,allRpts.expId]=unique(regexp(allRpts.name,'.*(?=_)','match','once'));
        nRpts=accumarray(allRpts.expId,1);
        if exist(h.meanFolder.String,'dir')>0
            meanFiles=struct2table(dir(h.meanFolder.String));
            meanExp=extractBefore(meanFiles.name,'.mat');
            existMeans=ismember(exps,meanExp);
            meanType={';yes)';';no)'};
            h.expList.String=strcat(exps,' (',num2str(nRpts),meanType(2-existMeans));
        else
            h.expList.String=strcat(exps,' (',num2str(nRpts),';nA)');
        end
        h.expList.Min=0;
        h.expList.Max=numel(h.expList.String);
        if numel(selExp)>0
            h.expList.Value=find(ismember(regexp(h.expList.String,'^.*(?= \(\d+;)','match','once'),selExp));
        else
            h.expList.Value=[];
        end
    end
    h.allRpts=allRpts;
    guidata(h.fig,h);
end

function calcMeans(a,b,type)
    h=guidata(a.Parent);
    GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
    load('/home/labs/barkailab/tamarge/Master/mat Files/promoterLengthsORF.mat','promoterLengthsORF');
    load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Divya/promType.mat','promType')
    if exist([h.meanFolder.String ,'/goodRepeats.xlsx'],'file')
        goodTable=readtable([h.meanFolder.String ,'/goodRepeats.xlsx']);
        if strcmp(type,'add')
            existFiles=dir([h.meanFolder.String ,'/*.mat']);
            selExp=find(~ismember(goodTable.Experiment,extractBefore({existFiles.name},'.mat')));
        elseif strcmp(type,'all')
            selExp=1:size(goodTable,1);
        elseif strcmp(type,'selection')
            selExp=find(ismember(goodTable.Experiment,extractBefore(h.expList.String(h.expList.Value),' (')));
        end
        for e=selExp(:)'
            selRpts=find(ismember(h.allRpts.name,strcat(goodTable.Experiment{e},'_',strsplit(goodTable.Repeat{e},';')','.mat')));
            c=0
            meanProfile=nan(sum(GP.chr_len),numel(selRpts));
            for r=selRpts'
                c=c+1;
                temp=matfile([h.allRpts.folder{r} '/' h.allRpts.name{r}]);
                meanProfile(:,c)=temp.normProfile;
                fn=fieldnames(temp);
                spName=fn(contains(fn,'sumprom','IgnoreCase',true));
                sumPromR(:,c)=temp.(spName{1});
            end
            meanProfile=mean(meanProfile,2,'omitnan');
            metaPro=metaProfilePromLenDivya(meanProfile,'promEnd','position','afterTss',0,'promLen',promoterLengthsORF);
            sumProm=squeeze(sum(metaPro,2,'omitnan'));
            sumProm(isnan(promoterLengthsORF),:)=NaN;
            clear temp
            temp.normProfile=meanProfile;
            temp.name=goodTable.Experiment{e};
            temp.sumProm=sumProm;
            temp.repeats=h.allRpts.name(selRpts);
            temp.name
            if numel(selRpts)>1
                temp.rptCorr=median(1-squareform(1-corr(sumPromR(promType<3,:),'rows','pairwise')));
            else
                temp.rptCorr=nan;
            end
            save([h.meanFolder.String ,'/', temp.name '.mat' ],'-struct','temp','-v7.3')
        end
    else
        action = 0
    end
end


%% these functions are for the check repeat figure

function checkRepeats(a,b)
    h=guidata(a.Parent);
    check.allRpts=h.allRpts;
    check.newExps=table (regexp(h.expList.String(h.expList.Value),'^.*(?= \(\d+;)','match','once'),repmat({''},numel(h.expList.Value),1),'VariableNames',{'Experiment','Repeat'})    
    check.fig=figure;
    check.crMap=subplot(1,3,1);
    check.scPlot=subplot(1,3,2);
    check.goodRpt=uicontrol('Style','listbox','String',{},'Units','normalized','Position',[0.6916 0.25 0.21 0.60],'Min',0,'Max',0,'Value',[],'CallBack',@updateGoodRpts);
    check.goodRptTitle=uicontrol('Style','text','String','Select good repeats','Units','normalized','Position',[0.6916 0.87 0.21 0.05]);
    check.prevExp=uicontrol('Style','pushbutton','String','Previous','Units','normalized','Position',[0.6916 0.12 0.05 0.05],'Callback',@(a,b)updateExp(a,b,-1));           
    check.nextExp=uicontrol('Style','pushbutton','String','Next','Units','normalized','Position',[0.75 0.12 0.05 0.05],'Callback',@(a,b)updateExp(a,b,1));
    check.saveRpt=uicontrol('Style','pushbutton','String','Save Good Repeats','Units','normalized','Position',[0.81 0.12 0.09 0.05],'Callback',@saveGoodRpt);   
    check.mainFig=h.fig;
    if size(check.newExps,1)>0
        check.expId=1;
        guidata(check.fig,check)
        refreshCheckFigure(check)
    else
        check.expId=[]
        guidata(check.fig,check)
    end       
end

function updateExp(a,b,dir)
    check=guidata(a.Parent);
    check.expId=check.expId+dir;
    if check.expId==0
        check.expId=size(check.newExps,1);
    elseif check.expId>size(check.newExps,1)
        check.expId=1;
    end
    guidata(check.fig,check)
    refreshCheckFigure(check)
end

function refreshCheckFigure(check)
    check.selRpts=find(ismember(regexp(check.allRpts.name,'.*(?=_)','match','once'),check.newExps.Experiment(check.expId)));
    c=0;
    check.sumProm=nan(6701,numel(check.selRpts));
    for r=check.selRpts'
        c=c+1;
        [check.allRpts.folder{r} '/' check.allRpts.name{r}]
        temp=matfile([check.allRpts.folder{r} '/' check.allRpts.name{r}]);
        tempFN=fieldnames(temp);
        spFields=tempFN(contains(tempFN,'sum','IgnoreCase',true)&contains(tempFN,'prom','IgnoreCase',true));
        if numel(spFields)>1
            spFields=spFields(contains(spFields,'new','IgnoreCase',true));
        end
        check.sumProm(:,c)=temp.(spFields{1});
        clear temp tempFN
    end      
    check.goodRpt.Min=0;
    check.goodRpt.Max=max(numel(check.selRpts),2); 
    check.goodRpt.String=strrep(extractBefore(check.allRpts.name(check.selRpts),'.mat'),'_',' ');  
    check.goodRpt.Value=find(ismember(regexp(check.allRpts.name(check.selRpts),'(?<=_)[A-Za-z0-9]+(?=.mat)','match','once'),strsplit(check.newExps.Repeat{check.expId},';')));    
    check.crIm=imagesc(check.crMap,corr(check.sumProm,'rows','pairwise'),[0.7 1]);
    colorbar(check.crMap)
    set(check.crMap,'XTick',1:numel(check.selRpts),'YTick',1:numel(check.selRpts),'YTickLabel',check.goodRpt.String,'XTickLabel',check.goodRpt.String,'XTickLabelRotation',25)
    title(check.crMap,sprintf('%s (%d/%d)',strrep(check.newExps.Experiment{check.expId},'_',' '),check.expId,size(check.newExps,1)))
    check.crIm.ButtonDownFcn=@(a,b)plotCorrelation(a,b,check);
    guidata(check.fig,check)
end

function plotCorrelation(a,b,check)
   xId=round(b.IntersectionPoint(1));
   yId=round(b.IntersectionPoint(2));
   scatter(check.scPlot,check.sumProm(:,xId),check.sumProm(:,yId),'.')
   xlabel(check.scPlot,check.goodRpt.String{xId})
   ylabel(check.scPlot,check.goodRpt.String{yId})
   title(check.scPlot,sprintf('%.2f',corr(check.sumProm(:,xId),check.sumProm(:,yId),'rows','pairwise')))
end

function updateGoodRpts(a,b)
    check=guidata(a.Parent);
    check.newExps.Repeat{check.expId}=strjoin(regexp(check.allRpts.name(check.selRpts(check.goodRpt.Value)),'(?<=_)[A-Za-z0-9]+(?=.mat)','match','once'),';');
    guidata(check.fig,check)
end

function saveGoodRpt(a,b)
    check=guidata(a.Parent);
    h=guidata(check.mainFig);
    if exist([h.meanFolder.String ,'/goodRepeats.xlsx'],'file')
        goodTable=readtable([h.meanFolder.String ,'/goodRepeats.xlsx']);
        [~,idx]=ismember(check.newExps.Experiment,goodTable.Experiment);
        goodTable.Repeat(idx(idx>0))=check.newExps.Repeat(idx>0);
        goodTable(end+1:end+sum(idx==0),:)=check.newExps(idx==0,:);
        action=1;
    else
        goodTable=check.newExps;
        action=2;
    end
    actionList={'updated','created'}
    writetable(goodTable,[h.meanFolder.String ,'/goodRepeats.xlsx']);
    msgbox(sprintf('%s %s',[h.meanFolder.String ,'/goodRepeats.xlsx'],actionList{action}))
end

    
