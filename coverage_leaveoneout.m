% CROSS VALIDATION SUCCESS RATE

% find % of points that are within the envelope (coverage interval) 
%inside has dimension 80 years of predictions x 79 models
clearvars
plotFlag=1;
%load model_full_priorsb_22_10_17.mat
load model_predictions_22_10_19.mat
load A20.mat
glm=runInfo.glm;
t2 = (1:122)';
inside=false(80,79);
validInt=43:122;
ct = 1;
for ct_mod = 1:size(A20,2)
    for ct_em = 1:size(A20(ct_mod).X,2)
        y(:,ct) = A20(ct_mod).X(1:42,ct_em);
        y2(:,ct) = A20(ct_mod).X(:,ct_em);
        ct = ct+1;
    end
end
load model_full_posteriors_22_10_14.mat
for run=1:79
    yhat=zeros(size(chainAll,1),length(t2));
    for ii = 1:size(chainAll,1)
        yhat(ii,:) = glm(chainAll(ii,:,run),t2);
    end
    pd1=fitdist(p1_priors(5000:end,run),'Normal');
    pd2=fitdist(p2_priors(5000:end,run),'Normal');
    pd3=fitdist(p3_priors(5000:end,run),'Normal');
    mlp(run,:)=[pd1.mean pd2.mean pd3.mean];
    
    y_full(:,run) = runInfo.glm(mlp(run,:),t2);
    y_full2(:,run)=smoothPH(y2(:,run),20);
    aCt=1;
    for alpha=0
        yint=prctile(yhat,[7 93]);
        yint_all(run,aCt,:,:)=yint;
        inside=zeros(1,80);
        vInt=validInt(y_full(validInt,run)>0.25);
        inside=zeros(1,length(vInt));
        inside(y_full(vInt,run)>yint(1,vInt)' & y_full(vInt,run)<yint(2,vInt)' & y_full(vInt,run)>0.25)=true;
        inPct(aCt,run)=sum(inside)./length(vInt);
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