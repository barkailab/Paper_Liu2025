result = readtable("msn2denovo.xlsx");
result.Sequence = nt2aa(result.Sequence);
fastawrite('msn2denovo_AA.fasta',table2struct(result))
msn2denovo = getMeta('disorder_scores_msn2.xlsx');

%% 
msn2denovo.DISO = cellfun(@(x) sum(x>0.5), msn2denovo.MetaPredict);
msn2denovo.DISOF = cellfun(@(x) sum(x>0.5)./numel(x), msn2denovo.MetaPredict);
msn2denovo.MeanPred = cellfun(@mean, msn2denovo.MetaPredict);

msn2denovo.Sample = cellfun(@(x) strrep(x,' ','_'),msn2denovo.Sample,'UniformOutput',false);
msn2denovo.Sample = cellfun(@(x) strrep(x,'_aN110_aQ61_aT55_aS74_','_aS74_aT55_aN110_aQ61_'),msn2denovo.Sample,'UniformOutput',false);

toremove = find(contains(msn2denovo.Sample,{'NQ','ST'}));
msn2denovo(toremove,:) =[];
%% strainid
for i= 1:height(msn2denovo)
    try
    msn2denovo.strainid(i) = min(find(ismember(strains.strain,msn2denovo.Sample{i})));
    catch
        msn2denovo.strainid(i) = FirstBatch.strainid(find(ismember(FirstBatch.Var1, msn2denovo.Sample{i}))); 
    end
end
msn2denovo.maxcorr=strains.max_corr(msn2denovo.strainid);
