function [out,sc4red]=mer_occupancy(profile,nmer,varargin)
ip=inputParser;
ip.addParameter('region','promoter');
ip.addParameter('window',20);
ip.addParameter('convert',false);
ip.addParameter('intGenes',1:6701);
ip.addParameter('badGenes',[]);
ip.addParameter('par',false);
ip.addParameter('method','normal');
ip.parse(varargin{:})
if ~ismember(nmer,[7])
    load('SC_genome.mat');
    clear sc4base
    ntBases={'A','C','G','T'};
    for i=1:17
        [~,sc4base{i}]=ismember(upper(SC_genome(i).Sequence'),ntBases);
        sc4base{i}=sc4base{i}-1;
    end
    clear SC_genome
    sc4base=cat(1,sc4base{:});
    kernel=4.^[0:nmer-1];
    %% create mers table
    motifVal=cellfun(@str2num,mat2cell(dec2base([0:4^nmer-1]',4,nmer),ones(4^nmer,1),ones(nmer,1)));
    nmerTable=table(mat2cell(cell2mat(ntBases(motifVal+1)),ones(4^nmer,1),nmer),[0:4^nmer-1]'+1,'VariableNames',{'seq','value'});
    nmerTable.rcValue=sum((3-motifVal).*kernel,2)+1;
    nmerTable.rcSeq=nmerTable.seq(nmerTable.rcValue);
    sc4nmer=conv(sc4base,kernel','same')+1;
    if ip.Results.convert & ismember(nmer,[5,7])
        load('motifs4ever.mat')
        if nmer==5
            nmerRed=nmerTable(ismember(nmerTable.seq,motifs4ever.motifSeqs5mer(motifs4ever.keptSeqs5mer)),:);
        else
            nmerRed=nmerTable(ismember(nmerTable.seq,motifs4ever.motifSeqs7mer(motifs4ever.keptSeqs7mer)),:);
        end
    else
        nmerRed=nmerTable(nmerTable.value<nmerTable.rcValue,:);
    end
    %sc4red=changem(sc4nmer,nmerRed.value,nmerRed.rcValue);
    %sc4red=changem(sc4red,1:size(nmerRed,1),nmerRed.value);
    sc4red=changem(sc4nmer,repmat([1:size(nmerRed,1)],1,2),[nmerRed.value;nmerRed.rcValue]');
    nmerRed.n=accumarray(sc4red,1,[size(nmerRed,1) 1],@sum,0);
    clear sc4nmer nmerTable sc4bases motifVal sc4base i ntBases kernel
else
    load(sprintf('genome%dmer.mat',nmer))
end
load('GP.mat');
chr_idx=[0 cumsum(GP.chr_len)];
promLen=700;
GP.gene_infoR64.status(cellfun('isempty',GP.gene_infoR64.status))={'Dubious'};
tssData='felixTss';
if ischar(ip.Results.region)
    if strcmpi(ip.Results.region,'promoter')
        gb=false(sum(GP.chr_len),1);
        for i=setdiff(ip.Results.intGenes,ip.Results.badGenes)
            if all(~isnan(GP.gene_infoR64.(tssData)(i,:)))
                if diff(GP.gene_infoR64.(tssData)(i,2:3))>0
                    prom=chr_idx(GP.gene_infoR64.(tssData)(i,1))+[GP.gene_infoR64.(tssData)(i,2)-promLen:GP.gene_infoR64.(tssData)(i,2)];
                else
                    prom=chr_idx(GP.gene_infoR64.(tssData)(i,1))+[GP.gene_infoR64.(tssData)(i,2):GP.gene_infoR64.(tssData)(i,2)+promLen];
                end
                gb(prom)=true;
            end
        end
        clear prom
        notGb=false(size(gb));
        for i=find(~contains(GP.gene_infoR64.status,'Dubious')| all(~isnan(GP.gene_infoR64.(tssData)),2))'
            if all(~isnan(GP.gene_infoR64.(tssData)(i,:)))
                cds=chr_idx(GP.gene_infoR64.(tssData)(i,1))+[min(GP.gene_infoR64.(tssData)(i,2:3)):max(GP.gene_infoR64.(tssData)(i,2:3))];
            else
                cds=chr_idx(GP.gene_infoR64.position(i,1))+[min(GP.gene_infoR64.position(i,2:3)):max(GP.gene_infoR64.position(i,2:3))];
            end
            notGb(cds)=true;
        end
        gb(notGb)=false;
    else
        gb=true(sum(GP.chr_len),1);
    end
else
    if numel(ip.Results.region)==sum(GP.chr_len)
        gb=acol(ip.Results.region);
    else
        gb=false(sum(GP.chr_len),1);
        gb(ip.Results.region)=true;
    end
end

nmerRed.n=accumarray(sc4red(gb),1,[size(nmerRed,1) 1],@sum,0);
%occ=movmean(profile,2*ip.Results.window+1,1);
score=zeros(size(nmerRed,1),size(profile,2));
score1=zeros(size(nmerRed,1),size(profile,2));
score2=zeros(size(nmerRed,1),size(profile,2));
if ip.Results.par
    parfor i=1:size(profile,2)
        if strcmp(ip.Results.method, 'normal')
            occ_i=movmean(profile(:,i),ip.Results.window,1);
        else
            %occ_i= max(movmean(profile(:,i),ip.Results.window,1) - movmean(profile(:,i),nmer,1),0);
            occ_i= max(movmean(profile(:,i),ip.Results.window,1) - movmean(profile(:,i),5,1).*5/ip.Results.window ,0);
        end
        score(:,i)=accumarray(sc4red(gb),occ_i(gb),[size(nmerRed,1) 1],@(x)mean(x,'omitnan'));
    end
else
    for i=1:size(profile,2)
        if strcmp(ip.Results.method, 'normal')
            occ_i=movmean(profile(:,i),ip.Results.window,1);
        else
            %occ_i= max(movmean(profile(:,i),ip.Results.window,1) - movmean(profile(:,i),nmer,1),0);
            occ_i= max(movmean(profile(:,i),ip.Results.window,1) - movmean(profile(:,i),5,1).*5/ip.Results.window ,0);
        end
        score(:,i)=accumarray(sc4red(gb),occ_i(gb),[size(nmerRed,1) 1],@(x)mean(x,'omitnan'));
    end
end
% score3=(score.*nmerRed.n-score2)./(nmerRed.n-1);
% for i=1:size(score,2)
%     subplot(4,4,i)
%     scatter(score(:,i),score1(:,i),[],log2(nmerRed.n),'.')
%     xlabel('mean')
%     ylabel('quantile')
%     subplot(4,4,i+4)
%     scatter(score(:,i),score2(:,i),[],log2(nmerRed.n),'.')
%      xlabel('mean')
%     ylabel('max')      
%     subplot(4,4,i+8)
%     scatter(score1(:,i),score2(:,i),[],log2(nmerRed.n),'.')
%     xlabel('90%')
%     ylabel('max')
%     subplot(4,4,i+12)
%     scatter(score(:,i),score3(:,i),[],log2(nmerRed.n),'.')
%     xlabel('mean')
%     ylabel('adj Mean')
% end
out.mers=nmerRed;
out.score=score;
end