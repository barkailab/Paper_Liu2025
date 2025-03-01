clear all
datafiles_rp=dir('/home/labs/barkailab/jingliu/Jingsamples/*.mat');
% 1 +bohdana strains
datafiles_BH=dir('/home/labs/barkailab/bohdana/001_sequencing_data/003_checProfiles');
%  get scramble and shift only
datafiles_rp = [datafiles_rp;datafiles_BH(find(contains({datafiles_BH.name},{'scramble','shift','MSN2_a','MSN2_v'})&~contains({datafiles_BH.name},'._')))];
%  get scramble and shift only

clearvars datafiles_BH

% 2 +divya strains
my_TF = unique(extractBefore({datafiles_rp.name},'_'));
datafiles_dk = dir('/home/labs/barkailab/divyakr/matlab/FinalForThesis/checProfiles');
%  for truncations: datafiles_dk(find(contains({datafiles_dk.name},'_d')&contains({datafiles_dk.name},my_TF)&~contains({datafiles_dk.name},{'_s','_r'})));
TFdescription = readtable('TFdescription.xlsx');
searchTerms = [TFdescription.FL;TFdescription.dkDBD; {'nonDBD'}];
searchTerms = searchTerms(~cellfun('isempty', searchTerms));

datafiles_rp = [datafiles_rp;datafiles_dk(find((contains({datafiles_dk.name}, ...
    searchTerms)&contains({datafiles_dk.name},my_TF))&~contains({datafiles_dk.name},{'_N','_s','_r','_p'})))];

msn2_files = dir("/home/labs/barkailab/jingliu/msn2_samples/*.mat");

% datafiles_rp = [datafiles_rp;msn2_files];

clearvars -except datafiles_rp msn2_files
%% if you want individual repeats
for i=1:numel(datafiles_rp)
    temp=load([datafiles_rp(i).folder '/' datafiles_rp(i).name]);
%        temp = load([datafiles_rp{i} '/' datafiles_rp(i).name]);
       if isfield(temp,'sumProm') && isfield(temp,'totalReads')
           JingsumProm(:,i)=temp.sumProm;
           sample{i}=strrep(datafiles_rp(i).name,'.mat','');
            totalReads(i) = temp.totalReads;
       elseif isfield(temp,'sumPromNew') && isfield(temp,'meta')
           JingsumProm(:,i)=temp.sumPromNew;
           sample{i}=strrep(datafiles_rp(i).name,'.mat','');
           try
           totalReads(i) = temp.meta.pairseUsed;
           catch
               totalReads(i) = temp.meta.uniquePairs;
           end
       elseif isfield(temp,'sumPromNew')
           JingsumProm(:,i)=temp.sumPromNew;
           sample{i}=strrep(datafiles_rp(i).name,'.mat','');
       elseif isfield(temp,'sumProm')
           JingsumProm(:,i)=temp.sumProm;
           sample{i}=strrep(datafiles_rp(i).name,'.mat','');
       end
      
       if isfield(temp,"merScoreNew")
            merscore(:,i) = temp.merScoreNew;
        else
            [out,sc4red]=mer_occupancy(temp.normProfile,7,'window',31);

            temp.merScoreNew = out.score;
            save([datafiles_rp(i).folder,'/',datafiles_rp(i).name],'-v7.3','-struct','temp')
            merscore(:,i) = temp.merScoreNew;
            temp=[];
        end
         
end


%% make msn2_samples -- attach to 
% to get files
msn2_want = {'MSN2_LIVFWY_to_Q.mat';'MSN2_DE_to_N.mat';'MSN2_FWY_toL.mat';'MSN2_KR_to_A.mat';'MSN2_N_to_G.mat';'MSN2_FWY_shift.mat'};
% to add in strains list
%msn2_nickname = {'MSN2_neg_toNQ'; 'MSN2_shift';'MSN2_FWY_toL';'MSN2_KR_toNQ';'MSN2_hyd_toNQ';'MSN2_NQ_toG'};
denovolist_new = {'aS74_aT55_aN110_aQ61_aD16_aE14_aF67',...
    'aS74_aT55_aN110_aQ61_aD16_aE14_aY67',...
    'aS74_aT55_aN110_aQ61_aD16_aE14_aL67',...
    'aS74_aT55_aN110_aQ61_aD16_aE14_aL130',...
    'aS114_aG89_aN43_aQ53_aD14_aE15_aF63',...
    'aS95_aG76_aN60_aQ68_aF63'};
denovolist_old = {'v2_wiF_wiDE','v2_wiY_wiDE','v2_wiL_wiDE','aST129_aNQ171_aED30_aL130','aS114_aNQ96_aG89_aED29_aF63','aS95_aNQ128_aG76_aF63'};
% msn2_want = [msn2_want; denovolist_new'];
% %% change bohdana sample name to be mine
% % newlist = [denovolist_new,'aS74_aT55_aN110_aQ61_aL67','aS74_aT55_aN110_aQ61_aD16_aE14'];
% % newoldlist = [denovolist_old,'v2_wiL_woDE','v2_woLFY_wiDE'];
% for i = find(ismember(regexp(sample, '^(.*?)(?=_\d)', 'match', 'once'),FirstBatch.Var1))
%      temp =find(ismember(FirstBatch.Var1, regexp(sample{i}, '^(.*?)(?=_\d)', 'match', 'once')));
%      sample{i}=strrep(sample{i},  FirstBatch.Var1{temp}, FirstBatch.newname{temp});
% 
% end
% General names
FirstBatch = load('FirstBatch');
newsample = sample;
for i = 1:numel(sample)
    id = find(ismember(FirstBatch.Var1,sample{i}));
    if ~isempty(id)
        sample{i} = FirstBatch.newname{id};
    end
end
 clearvars id
%% make sample table
% get TF list
TF=extractBefore(sample,'_');
% get sample table for repeats
sTable=table(sample',TF','VariableNames',{'strain_rp','TF'});

[a,b]=regexp(sTable.strain_rp,'_(?=\d+$)','once','split');
sTable = [sTable,cell2table(cat(1,a{:,1}),'VariableNames',{'strain','repeat_id'})];

% non necessary # of reads
sTable.reads = totalReads';
sTable.good = sTable.reads>1e6;
sTable.person = extractBefore(extractAfter({datafiles_rp.folder},'barkailab/'),'/')';
[strains,~,sTable.strainId]=unique(sTable.strain);

strains = strrep(strains,'_d2_','_d02_');
strains = strrep(strains,'700_d02_','700_d2_');
strains = strrep(strains,'_d95_','_d095_');
% add selected grammar and denovo msn2 
msn2_addstrains = msn2_files((contains({msn2_files.name},[msn2_want;{'_a'};{'_v'}])),:);

clearvars a b temp 
%% make sp for strains
datafiles_m = dir('/home/labs/barkailab/jingliu/Jingmeans/*.mat');
datafiles_m = [datafiles_m;msn2_addstrains];
load("TF_reference_7mer.mat")
load("TF_reference.mat")
FirstBatch = readtable('FirstBatch.xlsx');
%% get correlation mean and max for each strain 
load("Protype.mat")
strains = cell2table(strains,'VariableNames',{'strain'});
strains.mean_corr = zeros(height(strains),1);
strains.max_corr = zeros(height(strains),1);
strains.num_rp = zeros(height(strains),1);
% strains.bestrp = [];
% strains.bestrpP =[];

for i=1:height(strains)
    repeats=find(sTable.strainId==i);
    all_corr = [];
    if numel(repeats)>1
        pairs = nchoosek(repeats,2);
        all_corr = [];
        for p = 1:height(pairs)
            all_corr(1,p) = corr(JingsumProm(Protype<3,pairs(p,1)),JingsumProm(Protype<3,pairs(p,2)));
            all_corr(2,p) = corr(merscore(:,pairs(p,1)),merscore(:,pairs(p,2)));
        end
    %     if ~isempty(all_corr == 1)
    %         for p = find(==1)'
    %         sTable.delete(pairs(p,1)) = 'T';
    %         end
    %     end
        strains.pairs{i} = pairs;
        strains.all_corr{i} = all_corr;
        strains.mean_corr(i) = mean(all_corr(1,:));
        [strains.max_corr(i),tempid] = max(all_corr(1,:));
        strains.bestrps{i} = pairs(tempid,:);
        strains.people{i} = sTable.person(repeats);
        strains.numRps(i) = numel(repeats);
        if max(all_corr(1,:))<0.9 & isempty(find(ismember(strrep({datafiles_m.name},'.mat',''),strains.strain{i})))
            [maxcorr, id] = max(all_corr(1,:));
            [~,choos] = max(sTable.reads(pairs(id,:)));
            strains.bestrp(i) = pairs(id,choos); 
            strains.bestrpP{i} = sTable.person{pairs(id,choos)};
        end
    end
    strains.num_rp(i) = numel(repeats);
    if ~isempty(all_corr)
        strains.mean_corr_m(i) = mean(all_corr(2,:));
        strains.max_corr_m(i) = max(all_corr(2,:));
    end
end
clearvars all_corr pairs repeats
%% store best repeats for scramble_shift

tempss = find(contains(strains.strain,'scramble_shift'));
for i = tempss'
    temp = find(contains(strains.strain,{'scramble','shift'})&~contains(strains.strain,'scramble_shift') ...
        &contains(strains.strain,strains.TF{i}));
    strains.bestrps{i} = cat(1,strains.bestrps{temp});
end

%% add strains_m
for i = 1:height(strains)
     if any(find(ismember([FirstBatch.Var1,FirstBatch.newname],strains.strain{i})))
         curname = [strains.strain(i);FirstBatch.Var1{ismember(FirstBatch.newname,strains.strain{i})};FirstBatch.newname{ismember(FirstBatch.Var1,strains.strain{i})}];
        mid = find(ismember(strrep({datafiles_m.name},'.mat',''),curname));
     else
         mid = find(ismember(strrep({datafiles_m.name},'.mat',''),strains.strain{i}));
     end
    
    if ~isempty(mid)
        if numel(mid)>1
            mid = mid(~contains({datafiles_m(mid).folder},'msn2_samples'));
        end
        temp=load([datafiles_m(mid).folder '/' datafiles_m((mid)).name]);
        StrainSumProm(:,i) = temp.sumProm;
        
        if isfield(temp,'merScoreNew')
            Strain7mer(:,i) = temp.merScoreNew;
        else
            [out,sc4red]=mer_occupancy(temp.normProfile,7,'window',31);
        
            Strain7mer(:,i) = out.score;
            temp.merScoreNew = out.score;
            save([datafiles_m(mid).folder '/' datafiles_m(mid).name],'-v7.3','-struct','temp')
            temp=[];
        end
         if strains.max_corr(i)<0.01&&isfield(temp,'rptCorr')
             
             strains.max_corr(i) = temp.rptCorr;
         end
         strains.person{i} = extractBefore(extractAfter(datafiles_rp(mid).folder,'barkailab/'),'/');
    elseif isfield(TF_reference,strains.strain{i})
        StrainSumProm(:,i) = TF_reference.(strains.strain{i});
        Strain7mer(:,i) = TF_reference_7merscore.(strains.strain{i});
        strains.person{i} ='divyakr';
%     else
%         mid = find(ismember(strrep({msn2_addstrains.name},'.mat',''),strains.strain{i}));
%         temp=load([msn2_addstrains(mid).folder '/' msn2_addstrains(mid).name]);
%         StrainSumProm(:,i) = temp.sumProm;
%         Strain7mer(:,i) = temp.merScoreNew;
%         if isfield(temp,'rptCorr')
%             strains.max_corr(i) = temp.rptCorr;
%         end
    
    elseif find(ismember(sTable.strain,strains.strain{i})) & (strains.bestrp(i)~=0)
            StrainSumProm(:,i) = JingsumProm(:,strains.bestrp(i));
            Strain7mer(:,i) = merscore(:,strains.bestrp(i));
            strains.person{i} = extractBefore(extractAfter(datafiles_rp(strains.bestrp(i)).folder,'barkailab/'),'/');
    end

end
%% felix rpt corr
for i = 1:height(msn2_addstrains)
         mid = find(ismember(strains.strain,strrep({msn2_addstrains(i).name},'.mat','')));
     if ~isempty(mid)
         if strains.max_corr(mid)<0.01
         temp=load([msn2_addstrains(i).folder '/' msn2_addstrains(i).name]);
         strains.max_corr(mid) = temp.rptCorr;
         end
     end
end
%% change bohdana names
for i = find(ismember(strains.strain,FirstBatch.Var1))'
     temp =find(ismember(FirstBatch.Var1, strains.strain{i}));
     strains.strain{i}=strrep(strains.strain{i},  FirstBatch.Var1{temp}, FirstBatch.newname{temp});

end
%% add not in strains 
toadd = find(~ismember(strrep({datafiles_m.name},'.mat',''),strains.strain));
toadd = toadd(~(ismember(strrep({datafiles_m(toadd).name},'.mat',''),[FirstBatch.Var1,FirstBatch.newname])));
cur_strains = width(StrainSumProm);
for i = 1:numel(toadd)
    mid = toadd(i);
    temp=load([datafiles_m(mid).folder '/' datafiles_m(mid).name]);
    StrainSumProm(:,i+cur_strains) = temp.sumProm;
    Strain7mer(:,i+cur_strains) = temp.merScoreNew;
  strains.strain{i+cur_strains} = strrep(datafiles_m(mid).name,'.mat','');
  if isfield(temp,'rptCorr')
      strains.max_corr(i+cur_strains)=temp.rptCorr;
  end
end


%% add missing 7mers
checked = find(~any(Strain7mer));
DK_normP = dir('/home/labs/barkailab/divyakr/matlab/FinalForThesis/newFiles/structs/medianNorm/*.mat');

for i = checked
    id = find(ismember({DK_normP.name},[strrep(strains.strain{i},'_d02_','_d2_'),'.mat']));
    if ~isempty(id)  
        temp = load([DK_normP(id).folder '/' DK_normP(id).name]);
         [out,sc4red]=mer_occupancy(temp.medianNorm,7,'window',31);
        Strain7mer(:,i) = out.score;
    end
end
%% add TF info
strains.TF = extractBefore(strains.strain,'_');
strains.TF(ismember(strains.TF,'')) = cellfun(@(x) x(1:4), strains.strain(ismember(strains.TF,'')), 'UniformOutput', false);
[TF,~,strains.TFId]=unique(strains.TF);
TFdescription = readtable('TFdescription.xlsx');

%% FL & DBD list

DBDid = find(contains(strains.strain,TFdescription.DBD));
[~,id] = sort(strains.TFId(DBDid));
DBDid = DBDid(id);
FLid = find(ismember(strains.strain,TFdescription.FL));
[~,id] = sort(strains.TFId(FLid));
FLid = FLid(id);

clearvars id
%% add FACS abundance
% facssave = strains(find(strains.Abs_log2median>0),:);
% writetable(facssave,'FACSdata.xlsx')
facssave = readtable('FACSdata.xlsx');
for i = 1:height(facssave)
    id = find(ismember(strains.strain,facssave.strain{i}));
    strains.Abs_log2median(id) = facssave.Abs_log2median(i);
    strains.d_log2median(id) = facssave.d_log2median(i);
end

%% assign types
denovolist_new = {'aS74_aT55_aN110_aQ61_aD16_aE14_aF67'
'aS74_aT55_aN110_aQ61_aD16_aE14_aY67'
'aS74_aT55_aN110_aQ61_aD16_aE14_aL67'
'aS74_aT55_aN110_aQ61_aD16_aE14_aL130'
'aS114_aG89_aN43_aQ53_aD14_aE15_aF63'
'aS95_aG76_aN60_aQ68_aF63'};
grammarlist = {'hyd_toNQ','neg_toNQ','LIV_toF','FWY_toL','KR_toNQ','NQ_toG','scramble','shift','scramble_shift'};
pertubationlist = {'FL','DBD'};
pertubationlist(13:21) = grammarlist;
pertubationlist(31:36)= denovolist_new;
pertubationlist(12) = {'nonDBD'};
strains.type = 999*ones(height(strains),1);
strains.type(ismember(strains.strain, TFdescription.FL)) = 0;
strains.type(ismember(strains.strain, TFdescription.DBD)) = 40;
strains.type(ismember(strains.strain, TFdescription.dkDBD)) = 40;
strains.type(contains(strains.strain, 'nonDBD')) = 12;
strains.type(contains(strains.strain, {'hyd','MSN2_LIVFWY_to_Q'})) = 13;
strains.type(contains(strains.strain, {'neg','DE_to_N'})) = 14;
strains.type(contains(strains.strain, 'toF')) = 15;
strains.type(contains(strains.strain, 'toL')) = 16;
strains.type(contains(strains.strain, 'KR')) = 17;
strains.type(contains(strains.strain, {'NQ_toG','MSN2_N_to_G'})) = 18;
strains.type(contains(strains.strain, 'scramble')&~contains(strains.strain,'shift')) = 19;
strains.type(~contains(strains.strain, 'scramble')&contains(strains.strain,'shift')) = 20;
strains.type(contains(strains.strain, '_scramble_shift')) = 21;
for i = 1:6
    strains.type(ismember(cellfun(@(x) extractAfter(x,'_'), strains.strain, 'UniformOutput', false), denovolist_new{i})) = 30+i;
end

strains.type(contains(strains.strain,'MSN2_a')&~contains(strains.strain,denovolist_new)) = 30;

%% clear vars
clearvars allref checked facssave grammarlist i id maxcorr mid msn2_addstrains msn2_files msn2_want out p pertubationlist sample sc4red temp tempid tempss toadd totalReads
%% add abs
for i = find(contains(strains.strain,'_scramble_shift'))'
    id = find(strains.TFId == strains.TFId(i)&strains.type == 19);
    id2 = find(strains.TFId == strains.TFId(i)&strains.type == 20);
    strains.Abs_log2median(i) =strains.Abs_log2median(id);
    strains.d_log2median(i) = strains.d_log2median(id);
    tempm=[];
    tempm = strains.max_corr([id,id2]);
    strains.max_corr(i) = mean(tempm(tempm>0));
    tempm = strains.max_corr_m([id,id2]);
    strains.max_corr_m(i) = mean(tempm(tempm>0));
end
% strains.strain{height(strains)+1} = 'STP1_scramble_shift';
% temp = load('/home/labs/barkailab/jingliu/Jingmeans/STP1_scramble_shift.mat');
% StrainSumProm=[StrainSumProm,temp.sumProm];
clearvars rps mid temp maxcorr id choos
%% make WT motifs binding scores 
for i = [1:8,10:20,22:26]    
% for i=[17,19,24] %no nonDBD only
 %   for i =13, wt = 122, dbdonly = 127 in Jingmeans
    BStable = BStables.(TF{i});
    bsSur = 70;
    sel = (BStable.type==1);
    clearvars motifMat
    normProfile = zeros(12157105,1);
    sel_strain = find(ismember(strains.TF,TF{i}));
    sel_strain = sel_strain(ismember(sel_strain,[acol(drawlist);denovoinfo.strainid;FirstBatch.strainid]));
    c=0;
    for j = sel_strain'
        c=c+1;
        tempname = strains.strain{j};
        tempid = min(find(ismember({DK_normP.name},[tempname,'.mat'])));
        if isempty(tempid)
            if strains.bestrp(j) == 0
                if find(ismember(FirstBatch.Var1,tempname))
                    tempname = FirstBatch.newname{find(ismember(FirstBatch.Var1,tempname))};
                elseif find(ismember(FirstBatch.newname,tempname))
                    tempname = FirstBatch.Var1{find(ismember(FirstBatch.newname,tempname))};
                elseif find(ismember(TFdescription.DBD,tempname))
                    tempname = TFdescription.dkDBD{find(ismember(TFdescription.DBD,tempname))};
                end
                tempid = min(find(ismember({DK_normP.name},[tempname,'.mat'])));
                temp = load([DK_normP(tempid).folder '/' DK_normP(tempid).name]);
                
            else
                tempid =  strains.bestrp(j);
                temp = load([datafiles_rp(tempid).folder '/' datafiles_rp(tempid).name]);
            end
        else
            temp = load([DK_normP(tempid).folder '/' DK_normP(tempid).name]);
        end
            
        if isfield(temp,'medianNorm')
            normProfile(1:numel(temp.medianNorm),c) = temp.medianNorm;
        elseif isfield(temp,'normProfile')
            normProfile(1:numel(temp.normProfile),c) = temp.normProfile;
        end
    end
        motifMat=reshape(normProfile(acol(round(BStable.pos(sel))+[-bsSur:bsSur].*(BStable.dir(sel))),:),sum(sel),2*bsSur+1,width(normProfile));    
        motifBind=squeeze(mean(motifMat,2));
        promBaseBind=mean(normProfile(promoterIDXvecF==1,:),1);
        bindScore=motifBind./promBaseBind;
        [~,idwt] = sort(bindScore(:,1),'descend');  
        mean_promoter = squeeze(mean(motifMat(:,35:105,:),[1,2]));
        strains.BindScore(sel_strain)=mean_promoter./promBaseBind';
  
        % draw in vitro motif binding plot
        mean_motifMat = squeeze(mean(motifMat,1));
        rescaled = mean_motifMat./max(mean_motifMat,[],1);
        
%         clf
%         imagesc(rescaled')
%         xticks(0:35:141)
%         xticklabels([-70:35:70])
%         yticks(1:10)
%         yticklabels(labels)
%         colorbar
%         title([TF{i} ' bind around motif'])
%         saveas(gca, fullfile([path, TF{i},'motif_bind_oldway.jpg']))
%         clearvars rescaled gnames labels 
end    
%% some statistics
%4 SS
temp.a = strains.all_corr(drawlist(1:25,19));
temp.b = strains.max_corr(drawlist(1:25,19));
%4 charge polar
tempcorr = acol(maxcorr(1:25,[14,17,18])>0.9);
tempans = acol(corrDBD(1:25,[14,17,18]));
tempans = tempans(tempcorr);
%% for fig 4567
allcolors = [
    0.1216, 0.4667, 0.7059;  % Blue
    1.0000, 0.4980, 0.0549;  % Orange
    0.1725, 0.6275, 0.1725;  % Green
    0.8392, 0.1529, 0.1569;  % Red
    0.5804, 0.4039, 0.7412;  % Purple
    0.5490, 0.3373, 0.2941;  % Brown
    0.8902, 0.4667, 0.7608;  % Pink
    0.4980, 0.4980, 0.4980;  % Gray
    0.7373, 0.7412, 0.1333;  % Olive
    0.0902, 0.7451, 0.8118;  % Cyan
    0.8, 0.6, 0.6;  % Dark red or maroon
    0.3176, 0.8588, 0.8588;  % Light Teal (changed from Light Orange)
    0.7922, 0.6980, 0.8392;  % Light Purple
    1.0000, 0.8000, 0.2000;  % Yellow
    0.6941, 0.3490, 0.1569;  % Dark Brown
    0.8471, 0.7490, 0.8471;  % Lavender
    175/255,210/255,176/255;  % F
    0.0000, 0.0000, 0.0000;  % Black
    173/255,208/255,175/255;  % F
    0.33, 0.47, 0.13;  % Darker Deep Blue (Darker version of 19)
   255/255, 220/255, 151/255;  % L
    0.5, 0.43, 0.29;  % Darker Moss Green (Darker version of 27)

];