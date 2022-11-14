% CROSS VALIDATION SUCCESS RATE

% find % of points that are within the envelope (coverage interval) 
%inside has dimension 80 years of predictions x 79 models
clearvars
plotFlag=1;
load model_full_priorsb_22_11_10.mat
load model_predictions_22_11_11.mat
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
%load model_full_priorsb_22_11_10.mat
for run=1:79
    yhat=glmtimeseries(glm,squeeze(chainAll(:,:,run)),t2);
    mlp(run,:)=[median(p1_priors(:,run)) median(p2_priors(:,run)) median(p3_priors(:,run)) median(p4_priors(:,run))]; 
    y_full(:,run) = glmtimeseries(glm,mlp(run,:),t2);
    y_ful2(:,run)=smoothPH(y2(:,run),20);
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
        hold on
        title(['Run ' num2str(mn)])
        pause
    end 
end