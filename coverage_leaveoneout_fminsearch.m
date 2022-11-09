% CROSS VALIDATION SUCCESS RATE

% find % of points that are within the envelope (coverage interval) 
%inside has dimension 80 years of predictions x 79 models
clearvars
plotFlag=1;
%load model_full_priorsb_22_10_17.mat
load model_predictionsb_22_10_18.mat
load A20.mat
glm=runInfo.glm;
t2 = (1:122)';
inside=false(80,79);
validInt=43:122;
sigObs=0.5934; %std from estimate_obs_std4
ct = 1;
for ct_mod = 1:size(A20,2)
    for ct_em = 1:size(A20(ct_mod).X,2)
        y(:,ct) = A20(ct_mod).X(1:42,ct_em);
        y2(:,ct) = A20(ct_mod).X(:,ct_em);
        ct = ct+1;
    end
end
%load model_full_posteriors_22_10_14.mat
mu1=mean(mean(chainAll(:,1,:)));
mu2=mean(mean(chainAll(:,2,:)));
mu3=mean(mean(chainAll(:,3,:)));
mu4=mean(mean(chainAll(:,4,:)));
for run=1:79
    yhat=zeros(size(chainAll,1),length(t2));
    for ii = 1:size(chainAll,1)
        yhat(ii,:) = glm(chainAll(ii,:,run),t2);
    end
    loglike=   @(b) sum(-1.*log(normpdf(y2(:,run)-glm(b,t2),0,sigObs)));
   [bfit,~] = fminsearch(loglike, [mu1 mu2 mu3 mu4]);
    
    y_full(:,run) = glm(bfit,t2);
    aCt=1;
    for alpha=0
        yint=prctile(yhat,[6 94]);
        yint_all(run,aCt,:,:)=yint;
        inside=zeros(1,80);
        inside(y_full(validInt,run)>yint(1,validInt)' & y_full(validInt,run)<yint(2,validInt)')=true;
        inPct(aCt,run)=sum(inside)./80;
        aCt=aCt+1;
    end
end
valPct=nanmean(inPct,2);

if plotFlag
    for mn=1:79
        clf;
        plot(squeeze(yint_all(mn,1,:,:))');
        hold on; 
        plot(y_full(:,mn));
        hold on;
        plot(y2(:,mn),'.')
        title(['Run ' num2str(mn)])
        pause
    end 
end