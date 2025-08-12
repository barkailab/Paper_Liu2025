function result = getMeta(filename)

    MetaPred = readtable(filename,'NumHeaderLines',1);

    for i = 1:height(MetaPred)
        for j = 1:numel(MetaPred.Var3{i})-1
            MetaPredDisO{i}(j) = MetaPred.(['Var',num2str(j+3)])(i);
        end
    end
    result = table();
    result.Sample = cellfun(@(x) x(1:end), MetaPred.Var2,'UniformOutput',false);
    result.Sequence = cellfun(@(x) x(1:end), MetaPred.Var3,'UniformOutput',false);
    result.MetaPredict = MetaPredDisO';
end
