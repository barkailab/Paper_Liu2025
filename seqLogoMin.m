function [finalFreq]= seqLogoMin(MotSet,meanSignal_sorted,varargin)
ip=inputParser;
ip.addParameter('showFigs',false)
ip.addParameter('rcLoop',true)
ip.parse(varargin{:});
rcLoop=ip.Results.rcLoop;

for i = 1:length(MotSet)
    motRank (i) = round((meanSignal_sorted(i) / sum(meanSignal_sorted)),2) * 100;
end

motRank = int64(motRank);

for i = 2:length(MotSet)
    [~,forAlign] = nwalign(MotSet{1},MotSet{i},'glocal',true);
    forAlignMatch = sum(forAlign(2,:)=='|');
    [~,revAlign] = nwalign(MotSet{1},seqrcomplement(MotSet{i}),'glocal',true);
    revAlignMatch = sum(revAlign(2,:)=='|');
    if revAlignMatch > forAlignMatch
        MotSet{i} = seqrcomplement(MotSet{i});
    end
end

RCMotSet = cellfun(@seqrcomplement, MotSet,'UniformOutput', false);
MotSetOriginal = MotSet;
if rcLoop
    for z = 1:2^(length(MotSet)-1)
        [score(z),finalFreq] = wmCheck(MotSetOriginal, z, motRank, RCMotSet);
    end
    [v, idx] = max(score);
else
    idx = 512;
end
[~,finalFreq] = wmCheck(MotSetOriginal, idx, motRank, RCMotSet);
if ip.Results.showFigs
    mySeqLogo(finalFreq);
end
end


function [score,finalFreq] = wmCheck(MotSetOriginal, z, motRank, RCMotSet)
    MotSet = MotSetOriginal;
    idxRev = regexp(dec2base(z-1,2,length(MotSet)-1),'0');
    idxRev = idxRev+1;
    MotSet(idxRev) = RCMotSet(idxRev);
    formatAlign2=[];
    c=0;
    for f = 1:length(motRank )
        for i = 1:motRank (f);
            c=c+1;
            formatAlign2(c).Header=num2str(c);
            formatAlign2(c).Sequence = MotSet{f};
        end
    end
%     j=0;
%     for i = 1:2:sum(motRank )*2-1
%         j=j+1;
%         finalMatforAlign{i} = ['>',num2str(j)];
%     end
%     e = 0;
%     for f = 1:length(motRank )
%         currSize = motRank (f);
%         for i = 1:currSize
%             finalMatforAlign{int64((e+i)*2)} = MotSet{f};
%         end
%         e = e + motRank (f);
%     end
    finalAlign = multialign(formatAlign2,'terminalGapAdjust',true);
    sumA = [];
    sumG = [];
    sumC = [];
    sumT = [];
    totAGCT = [];
    
    for h= 1:length(finalAlign(1).Sequence)
        currPos = char();
        for i = 1:length(finalAlign)
            currPos = [currPos,char(finalAlign(i).Sequence(h))];
            sumA(h) = sum(currPos=='A');
            sumC(h) = sum(currPos=='C');
            sumG(h) = sum(currPos=='G');
            sumT(h) = sum(currPos=='T');
            totAGCT(h) = length(regexp(currPos,'[AGCT]'));
        end
    end
    
    sumAfreq = sumA ./ totAGCT;
    sumGfreq = sumG ./ totAGCT;
    sumCfreq = sumC ./ totAGCT;
    sumTfreq = sumT ./ totAGCT;
    finalFreq = [sumAfreq;sumCfreq;sumGfreq;sumTfreq];
    
    %smsum = smooth(totAGCT);
    %Positions2Plot=find(smooth(totAGCT)>max(smooth(totAGCT))/1.8);
    motLen=numel(MotSet{1});
    if length(totAGCT) > motLen
    [~,iMax] = max(movsum(totAGCT,motLen,'Endpoints','shrink'));
    Positions2Plot = [iMax-(motLen-1)/2:iMax+(motLen-1)/2];
    [wm,~] = seqlogo(finalFreq(:,Positions2Plot(1):Positions2Plot(end)),'DisplayLogo',false);
    score = sum(sum(wm{2}.*totAGCT(Positions2Plot(1):Positions2Plot(end))));
    finalFreq = finalFreq(:,Positions2Plot(1):Positions2Plot(end));
    else
    [wm,~] = seqlogo(finalFreq,'DisplayLogo',false);
    score = sum(sum(wm{2}.*totAGCT));
        
%     if length(finalAlign(1).Sequence)>12
%         score = 0;
%     end
    end
end
