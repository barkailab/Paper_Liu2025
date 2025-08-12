%% get all PFMs
% DK_normP is combined files directory of data from this study and Kumar et al., 2023, Molecular
% Cell.
normProfile = [];
for i = 1:height(DK_normP)
        temp = load([DK_normP(i).folder,'/',DK_normP(i).name]);
        if isfield(temp,'temp')
            normProfile = temp.temp.normProfile;
        elseif isfield(temp,'medianNorm')
            normProfile = temp.medianNorm;
        else
            normProfile = temp.normProfile;
        end
        % calculate unbiased 7-mer motif appearance
        [out,sc4red]=mer_occupancy(normProfile,7,'window',21);
        [~,top10]=maxk(out.score,10);
        out.mers(top10,:);
        % calculate PFM
        [finalFreq]= seqLogoMin(out.mers.seq(top10),out.score(top10),'showFigs',false);
        PFM{i} = finalFreq;
        finalFreq=[];  
end 
save('PFM_mean.mat','PFM','DK_normP')
clearvars temp normProfile out sc4red top10
%% compare motifs
TFdescription=readtable('TFdescription.xlsx');
TFdescription.TFi = cellfun(@(x) [x(1) lower(x(2:end))], TFdescription.TF,'UniformOutput',false);
RefPWM = load_meme_pwm('Saccharomyces_cerevisiae_CISBP2.meme');
% load TF in vitro motif pwm
for i = 1:height(TFdescription)
    if ~isempty(TFdescription.CISBP_Id{i})
    pwmt{i} = RefPWM((ismember({RefPWM.name},TFdescription.CISBP_Id{i}))).matrix';
    end
end
clearvars RefPWM
%% get respective TFs of my samples
DK_TFs = cellfun(@(x) extractBefore(x,'_'),{DK_normP.name},'UniformOutput',false);
tempPl = find(cellfun(@isempty, DK_TFs));
DK_TFs(tempPl) = cellfun(@(x) strrep(x,'.mat',''),{DK_normP(tempPl).name},'UniformOutput',false);
clearvars tempPl
%% get corresponding strainid of these samples in the main 'strains.mat'
for i = 1:height(DK_normP)
    tempid = min(find(ismember(strains.strain,strrep(DK_normP(i).name,'.mat',''))));
    if ~isempty(tempid)
        DK_normP(i).strainid = tempid;
    else
        DK_normP(i).strainid = 0;
    end
end
clearvars tempid
%% compare PWM of each de novo motif and all in vitro motifs
ComPFMRef=zeros(height(DK_normP),numel(pwmt));
RefTFid=zeros(height(DK_normP),i);
for i=1:numel(DK_TFs)
    TFid = find(ismember(TFdescription.TF, DK_TFs{i}));
    if ~isempty(TFid)
        for j = 1:numel(pwmt)
            if ~isempty(pwmt{j})
                ComPFMRef(i,j) = comparePWM(PFM{i},pwmt{j});
            end
        end
        RefTFid(i) = TFid;
    end
end
TFid = [];

%% Collect motifs comparision between each sample and respective in-vitro motifs, and with all other in-vitro motifs
sel = [1:8,10:26]; % Hot1 has no in vitro motifs
withref = zeros(numel(DK_TFs),1);
withother = cell(numel(DK_TFs),1);
for i = 1:numel(RefTFid)
    withref(i) = ComPFMRef(i,sel(sel==RefTFid(i)));
    if ComPFMRef(i,sel(sel==RefTFid(i)))~=0
    % take out other TFs similar to the reference
    withother{i} = ComPFMRef(i,sel(sel~=RefTFid(i)));
    else
        withref(i)=[];
        withother{i}=[];
    end
end


%% Collect denovo motifs of mutants with respective FL; and with in-vitro
clearvars UnbMCI UnbMCJ
n=0;
for i = 1:numel(PFM)
    n=n+1;
    if DK_normP(i).strainid>0&~isempty(withref{i})
        UnbMCI(n) = withref{i}; %unbiased motifs compared to in-vitro
        wtid = find(ismember([DK_normP.strainid],drawlist(strains.TFId(DK_normP(i).strainid,1))));
        if ~isempty(wtid)
            if any(ismember(wtid,706))
                wtid = 706;
            end
            UnbMCJ(n) = comparePWM(PFM{i},PFM{wtid});%unbiased motifs compared to unbiased motif of respective FL
        end
    end
end
