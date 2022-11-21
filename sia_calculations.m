%Script to generate calculations for SIA manuscript
%
% Ted Amdur
% 11/21/22

clearvars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set to 1 to calculate ice-free dates given observations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iceFree=1;
if iceFree
    load obs_pred_22_11_16long.mat %Insert name of saved prediction one desires to use
    glm=obsInfo.glm;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set to 1 to calculate ice-free range from CMIP6 simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmip6=0;
if cmip6
    load mat_files/A20.mat
    load model_full_priors_22_11_14.mat
    glm=runInfo.glm;
    sigObs=0.6511; %std from estimate_obs_std4
    meanObs=0; %mean from estimate_obs_std4
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set to 1 to calculate parameter ranges of prior distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
priorStats=0;
if priorStats
    load model_full_priors_22_11_14.mat
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set to 1 to calculate parameter ranges of posterior distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
posteriorStats=0;
if posteriorStats
    load obs_pred_22_11_16long.mat
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Vectors used by multiple calculations:
yrObs=(1:42)'; %Years of observations to be used
t = yrObs;
tAll = (1:122)';
stYr=1978;

if iceFree
    yhat=glmtimeseries(glm,chain,tAll);
    %Find closest year for each draw
    [~,yr]=min(abs(yhat-1),[],2);
    yr=yr+stYr; %Change to year
    iceFreeDates=prctile(yr,[2.5 10 50 90 97.5]);
end
if cmip6
    ct = 1;
    for ct_mod = 1:size(A20,2)
        for ct_em = 1:size(A20(ct_mod).X,2)
            y2(:,ct) = A20(ct_mod).X(:,ct_em);
            ct = ct+1;
        end
    end
    %load model_full_posteriors_22_10_14.mat
    mu1=mean(p1_priors(:));
    mu2=mean(p2_priors(:));
    mu3=mean(p3_priors(:));
    mu4=mean(p4_priors(:));
    for ii=1:79
        loglike=  @(b) sum(-1.*log(normpdf(y2(:,ii)-glm(b,tAll),meanObs,sigObs)));
        [bfit,~] = fminsearch(loglike, [mu1 mu2 mu3 mu4]);
        yhat(ii,:)= glmtimeseries(glm,bfit,tAll);
    end
    [~,yrF]=min(abs(yhat-1),[],2);
    yrF=yrF+stYr; %Change to year
    cmipIceFree=prctile(yrF,[2.5 10 50 90 97.5]);
end
if priorStats
    A0=prctile(p1_priors(:),[2.5 10 50 90 97.5]);
    alpha=prctile(p2_priors(:),[2.5 10 50 90 97.5]);
    t0=prctile(p3_priors(:),[2.5 10 50 90 97.5]);
    nu=prctile(p4_priors(:),[2.5 10 50 90 97.5]);
    
end
if posteriorStats
    A0=prctile(chain(:,1),[2.5 10 50 90 97.5]);
    alpha=prctile(chain(:,2),[2.5 10 50 90 97.5]);
    t0=prctile(chain(:,3),[2.5 10 50 90 97.5]);
    nu=prctile(chain(:,4),[2.5 10 50 90 97.5]);
end
