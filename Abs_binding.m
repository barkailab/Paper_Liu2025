%% start with normProfile
function [motifMat,sumMat,TF]=absbinding(temp)
    normProfile = temp.normProfile;
    if size(normProfile,1)==1;
        normProfile= normProfile';
    end
    
    if size(normProfile,1)<12157105
        normProfile(size(normProfile,1)+1:12157105)=0;
    end
    load('BStables_1e3.mat')
    load('promoterIDXvecFimp.mat')
    TF = extractBefore(temp.name,'_');
    if ismember(TF,{'x','X'})|isempty(TF)
        TF = 'MSN2';
    end
    sumMat = zeros(4,1);
    if isfield(BStables,TF)
        BStable = BStables.(TF);
        bsSur = 30;
        sel = (BStable.type==1);
    
        motifMat=reshape(normProfile(acol(round(BStable.pos(sel))+[-bsSur:bsSur].*(BStable.dir(sel))),:),sum(sel),2*bsSur+1,width(normProfile));
    
        sumMat(1) = sum(motifMat(:));
    else
        motifMat=0;
    end
    sumMat(2) = sum(normProfile(promoterIDXvecF==1),1);%AllpromoterMat
    sumMat(3) = sum(normProfile(promoterIDXvecF==5),1);%mitochondira mat
    sumMat(4) = sum(normProfile(promoterIDXvecF==4),1);%TeloMat 
    
end