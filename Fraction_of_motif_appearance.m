%% get number of motifs in each promoter
load('BStables_1e3.mat')
load('promoterLengthsORF')
MotifNsumProm = zeros(6701,height(TFdescription));
for i= 1:height(TFdescription)
    if isfield(BStables,TFdescription.TF{i})
        BStable = BStables.(TFdescription.TF{i});
        sel = (BStable.type==1);
        tempNorm = zeros(12157105,1);
        tempNorm(round(BStable.pos(sel)))=1;
        metaPro= metaProfilePromLenDivya(tempNorm,'promEnd','position','afterTss',0,'promLen',promoterLengthsORF);
              sumProm=squeeze(sum(metaPro,2,'omitnan'));
        sumProm(isnan(promoterLengthsORF),:)=NaN;
        MotifNsumProm(:,i) = sumProm;
        sumProm = [];
    end
end
clearvars sel tempNorm metaPro BStables
%% Access all targets
load('JingSumProm.mat')
% log for msn2 & my pertubations
log_sum = log2(StrainSumProm+700); 
zscore= (log_sum-mean(log_sum,1,'omitnan'))./std(log_sum,[],1,'omitnan');
z_log = zscore>3; 
% interested strains
% TFId 27 is with absence of DBD, and TFId = 10 is Hot1, does not have invitro motifs  
IntStrains = [find(strains.TFId<27&strains.TFId>0&strains.TFId~=9&strains.max_corr>=0.9)]; 
% if youw ant de novo ismember(strains.type,[30:36])
motifYesMat=zeros(sum(Protype<3),width(z_log));
motifYesMat(:,IntStrains) = MotifNsumProm(Protype<3,strains.TFId(IntStrains))>0;

% sum number of z_log 
MotifTarFrac = sum(motifYesMat.*z_log(Protype<3,:))./sum(z_log(Protype<3,:));
sigmotifNocare = sum(motifYesMat)./sum(Protype<3);

figure
histogram(MotifTarFrac(IntStrains),'Normalization','probability')
hold on
histogram(sigmotifNocare(IntStrains),'Normalization','probability')
legend({'Targets(z>3)','All promoters(n=5358)'})
xlabel('Fraction of promoters with motifs')
ylabel('Probability')

%% full gradient changes
for i = 1:height(pertubationinfo)
    PerbMatMotifFra(i,1:numel(pertubationinfo.cursels{i})) = MotifTarFrac(denovoinfo.strainid(sel(pertubationinfo.cursels{i})));
end
