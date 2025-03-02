%% get all indivdual repeats
clear all
datafiles = dir('Y:\Jingsamples\*.mat'); % all files ending with mat

for i=1:numel(datafiles)
    % numel is # of elements
    temp=load([datafiles(i).folder '/' datafiles(i).name]);
    if isfield(temp,'sumProm')
        sumProm(:,i)=temp.sumProm;
        sample{i}=strrep(datafiles(i).name,'.mat','');
        totalReads(i) = temp.totalReads;
    elseif isfield(temp,'sumPromNew') & isfield(temp,'meta')
        sumProm(:,i)=temp.sumPromNew;
        sample{i}=strrep(datafiles(i).name,'.mat','');
        totalReads(i) = temp.meta.pairseUsed;
    else
        sumProm(:,i)=temp.sumPromNew;
        sample{i}=strrep(datafiles(i).name,'.mat','');
    end
    
end

[a,b]=regexp(sample,'_(?=[0-9]$)','once','split');
% separate name and repeat_ID
sampleTable=cell2table(cat(1,a{:}));

% stacking all strain names into a table
[strains_mean,~,sampleTable.strainId]=unique(sampleTable.Var1);
sampleTable.strain_rp = sample';

%% look at individual repeats
for i=1:numel(interestedstrains)
    repeats=find(ismember(sampleTable.Var1,interestedstrains{i}));
%       interestedstrainsall.allrp{i} = repeats;
%  end
    imagesc(corr(JingsumProm(:,repeats),'rows','pairwise'),[0.75 1]) % rows pairwise ignores NaNs
    set(gca,'Xtick',[1:numel(repeats)],'Ytick',[1:numel(repeats)],'Yticklabel',sampleTable.Var2(repeats),'xticklabel',string(repeats));
    title(strrep(interestedstrains{i},'_',' '))
    colorbar()
    pause()
end
% here save reproducible repeats in a excel called goodRepeats.xlsx
%%
load('promoterLengthsORF.mat')
load('Protype')
%% make means
goodRepExcel=readtable('goodRepeats.xlsx','ReadVariableNames',true);
exsitmean = dir('/home/labs/barkailab/jingliu/Jingmeans/*.mat');
for i=1:height(goodRepExcel)
    if isempty(find(ismember(strrep({exsitmean.name},'.mat',''),goodRepExcel.name{i})))
        selRepeats = strsplit(goodRepExcel.Repeat{i},";");
%         goodRepeats=find(ismember(sampleTable.Var1,goodRepExcel.Experiment{i})&ismember(sampleTable.Var2,selRepeats));
goodRepeats =str2double(selRepeats);
    if numel(goodRepeats)>1
            c=0;
            repeatProfile=[];
            for r=goodRepeats
                c=c+1;
                temp=load([datafiles_rp(r).folder '/' datafiles_rp(r).name]);
                repeatProfile(:,c)=temp.normProfile;
                repeatName{c}=strrep(datafiles_rp(r).name,'.mat','');
                repeatsp(:,c) = temp.sumProm;
            end
            meanProfile=mean(repeatProfile,2);
            [out,sc4red]=mer_occupancy(meanProfile,7,'window',31);

            
        metaPro=metaProfilePromLenDivya(meanProfile,'promEnd','position','afterTss',0,'promLen',promoterLengthsORF);
            sumProm=squeeze(sum(metaPro,2,'omitnan'));
            sumProm(isnan(promoterLengthsORF),:)=NaN;
            meanData.normProfile=meanProfile;
            meanData.sumProm=sumProm;
            meanData.repeats=repeatName;
            meanData.name=goodRepExcel.name{i};
            meanData.merScoreNew = out.score;
            meanData.rptCorr = corr(repeatsp(Protype<3,:),'rows','pairwise');
            save(['/home/labs/barkailab/jingliu/Jingmeans/',meanData.name ,'.mat'],'-v7.3','-struct','meanData')
            clear meanData sumProm metaPro c repeatProfile repeatName meanProfile temp
    end
    else
        goodRepExcel.name{i}
    end

end
