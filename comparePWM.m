function [disPwm, PWM1, PWM2, aF, alignCheck] = comparePWM(PWM1,PWM2,varargin)
ip = inputParser;
ip.addParameter('method', 'Euclidean');
ip.addParameter('coreLength', []);
ip.parse(varargin{:});

PWM1 = PWM1./sum(PWM1);
PWM2 = PWM2./sum(PWM2);

nPos1 = size(PWM1,2);
nPos2 = size(PWM2,2);

[~,maxIdx1] = max([PWM1; repmat(0.5,1,size(PWM1,2))]);
[~,maxIdx2] = max([PWM2; repmat(0.5,1,size(PWM2,2))]);

NT = {'A','C','G','T','N'};

cons1 = cat(2,NT{maxIdx1});
cons2 = cat(2,NT{maxIdx2});

aF = localalign(cons1, cons2, 'Alphabet', 'NT');
aR = localalign(cons1, seqrcomplement(cons2),'Alphabet', 'NT');

if max([aF.Score,0]) < aR.Score 
   PWM2 = fliplr(flipud(PWM2));
   cons2 = seqrcomplement(cons2);
   aF = aR;
end

clear aR
startdiff = diff(aF.Start);

if numel(startdiff)>0
    if startdiff > 0
        PWM1 = [repmat(0.25, 4, startdiff) , PWM1];
    elseif startdiff < 0
        PWM2 = [repmat(0.25, 4, abs(startdiff)) , PWM2];
    end
    
    PWM1 = [PWM1, repmat(0.25, 4, length(PWM2)-length(PWM1))];
    PWM2 = [PWM2, repmat(0.25, 4, length(PWM1)-length(PWM2))];
    
    disPwmVec = distanceCalculator(PWM1, PWM2, 'method',ip.Results.method);
    if numel(ip.Results.coreLength) == 0
        coreLength = max(5, length(aF.Alignment{1}));
    else
        coreLength = ip.Results.coreLength;
    end
    disPwm = min(movmean(disPwmVec,coreLength, 'Endpoints','fill'));
    alignCheck = diff(aF.Start) - diff(aF.Stop);
    
else
    disPwm = 1;
end

% 
% 
% if size(aF.Alignment{:},2)>=5 
%     PWM1 = PWM1(:,aF.Start(1):aF.Stop(1));
%     PWM2 = PWM2(:,aF.Start(2):aF.Stop(2));
% else
%     midAl1 = floor((aF.Start(1)+aF.Stop(1))/2);
%     takePos1 = zeros(size(PWM1(1,:)));
%     takePos1(midAl1)=1;
%     [~,takePos1]=maxk(smoothdata(takePos1,'gaussian',10),5)
%     PWM1 = PWM1(:,sort(takePos1));
%     
%     midAl2 = floor((aF.Start(2)+aF.Stop(2))/2);
%     takePos2 = zeros(size(PWM2(1,:)));
%     takePos2(midAl2)=1;
%     [~,takePos2]=maxk(smoothdata(takePos2,'gaussian',10),5)
%     PWM2 = PWM2(:,sort(takePos2));
%     
%     fprintf('Alignment length was shorter 5\n')
% 
% end
% if startdiff > 0
%     PWM1 = [repmat(0.25, 4, startdiff) , PWM1];
%     PWM2 = [PWM2, repmat(0.25, 4, startdiff)];
% elseif startdiff < 0
%     PWM2 = [repmat(0.25, 4, abs(startdiff)) , PWM2];
%     PWM1 = [PWM1, repmat(0.25, 4, abs(startdiff))];
%     
% end    
% disPwm = distanceCalculator(PWM1, PWM2, 'method',ip.Results.method);   
% alignCheck = diff(aF.Start) - diff(aF.Stop);

end


function disPWMVec = distanceCalculator(pwm1, pwm2, varargin)
ip = inputParser;
ip.addParameter('method', 'Euclidean');
ip.parse(varargin{:});
   nPos = size(pwm1, 2);

   if strcmp(ip.Results.method, 'Euclidean')
       disPWMVec = (1/(sqrt(2)))*(diag(pdist2(pwm1',pwm2')));
   elseif strcmp(ip.Results.method, 'Correlation')
       corrMat = corr(pwm1-0.25,pwm2-0.25);
       corrMat(isnan(corrMat)) = 0;
       disPWMVec = 1-(diag(corrMat));
%        disPWMVec = min(1, disPWMVec);
       
   end

end

function disPWM = distanceCalculatorOld(pwm1, pwm2, varargin)
ip = inputParser;
ip.addParameter('method', 'Euclidean');
ip.parse(varargin{:});
   nPos = size(pwm1, 2);

   if strcmp(ip.Results.method, 'Euclidean')
       
       disPWM = (1/(sqrt(2)*nPos))*sum(diag(pdist2(pwm1',pwm2')));
   elseif strcmp(ip.Results.method, 'Correlation')
       corrMat = corr(pwm1-0.25,pwm2-0.25);
       corrMat(isnan(corrMat)) = 0;
       disPWM = 1-(1/nPos)*sum(diag(corrMat));
       disPWM = min(1, disPWM);
       
   end

end