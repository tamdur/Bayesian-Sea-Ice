function [chainAll,y_full,runInfo]= posterior_estimate_cmip(saveStr)
%Create posterior estimates for general logistic function parameters using
%CMIP6 model runs of sea ice area
%
% Ted Amdur, Charlotte Dyvik Henke
% 10/13/22

%Define interval examined, slice sampler parameters
yrsSampled=(1:42)'; %Observations to be sampled 
nRuns=79; %Number of runs to use
leaveOneOut=1; %1 to leave out model run from prior dist, 0 otherwise
drawPriors=0; %1 to calculate mle's of full glm from each CMIP6 run, 0 otherwise
nsamples=10000;burn=1000;thin=5;%Number of draws; number of draws to first burn;thinning parameter
iF=1; %inflation factor
chainAll=zeros(nsamples,4,nRuns);
sigObs=0.5934.*iF; %std from estimate_obs_std4
meanObs=-0.2681; %mean from estimate_obs_std4

%Define the model
glm  = @(b,t) b(1).*(1-1./(1+exp(-b(2)*(t-b(3)))).^(1/b(4))); % SIA sigmoid

%Load CMIP6 model data
load A20.mat

% creating array of model values
t = yrsSampled;
tAll=(1:122)';
ct = 1;
y=zeros(length(t),nRuns);
y_full=zeros(length(tAll),nRuns);
for ct_mod = 1:size(A20,2)
    for ct_em = 1:size(A20(ct_mod).X,2)
        y(:,ct) = A20(ct_mod).X(t,ct_em);
        y_full(:,ct) = A20(ct_mod).X(tAll,ct_em);
        ct = ct+1;
    end
end

if drawPriors
    [p1_priors,p2_priors,p3_priors,~]= sigmoid_params(A20);
else
    load model_full_priorsb_22_10_17.mat
end

% storing most likely 1st round priors
y_hat_full=zeros(length(tAll),nRuns);
mlp=zeros(nRuns,4);

tElapsed=0;
if ~leaveOneOut
    p1=fitdist(p1_priors(:),'Lognormal');
    p2=fitdist(p2_priors(:),'Lognormal');
    p3=fitdist(p3_priors(:),'Normal');
    p4=fitdist(p4_priors(:),'Half Normal');
end

for ii = 1:nRuns
    tic
    %Create a prior that leaves out information from that model run
    incl=1:nRuns;
    if leaveOneOut
        incl=incl(incl~=ii)';
        p1p=p1_priors(:,incl);
        p2p=p2_priors(:,incl);
        p3p=p3_priors(:,incl);
        p4p=p4_priors(:,incl);
        p1=fitdist(p1p(:),'Lognormal');
        p2=fitdist(p2p(:),'Lognormal');
        p3=fitdist(p3p(:),'Normal');
        p4=fitdist(p4p(:),'Half Normal');       
    end
    mu1=p1.mean;mu2=p2.mean;mu3=p3.mean;mu4=p4.mean;
    sig1=p1.std;sig2=p2.std;sig3=p3.std;sig4=p4.std;
    
    %Get the fit over the full simulation for each CMIP run    
    pr1=fitdist(p1_priors(:,ii),'Lognormal');
    pr2=fitdist(p2_priors(:,ii),'Lognormal');
    pr3=fitdist(p3_priors(:,ii),'Normal');
    pr4=fitdist(p4_priors(:,ii),'Half Normal');
    mlp(ii,:)=[pr1.mean pr2.mean pr3.mean pr4.mean];
    y_hat_full(:,ii) = glm(mlp(ii,:),tAll);
    % ====================================================================== %
    % set initial values, mean, sigma.
    %b10 - A_0          %b20  - alpha        %b30 - t0
    init    = [mu1                mu2        mu3        mu4];
    %mu1                %mu2                 %mu3
    mu      = [mu1                mu2        mu3        mu4];
    %sigma1             %sigma2              %sigma3
    sigma   = [sig1.*iF           sig2.*iF   sig3.*iF   sig4.*iF];
    
    %Then, define the prior pdfs for the parameters
    prior1= @(b1) pdf(p1,b1);
    prior2= @(b2) pdf(p2,b2);
    prior3= @(b3) pdf(p3,b3);
    prior4= @(b4) pdf(p4,b4);
    
    %Write function for likelihood, posterior
    logpost=   @(b) sum(log(normpdf(y(:,ii)-glm(b,t),meanObs,sigObs)))+log(prior1(b(1)))+...
        log(prior2(b(2)))+log(prior3(b(3)))+log(prior4(b(4)));
    
    %Run slice sampler
    [chain,~]=slicesample(init,nsamples,'logpdf',logpost,'thin',thin,'burnin',burn);
    toc %Display elapsed time spent sampling
    
    % store values for each parameter
    chainAll(:,:,ii) = chain;
    nElapsed=toc;
    tElapsed=tElapsed+nElapsed;
    disp(['Run ' num2str(ii) ' posterior complete. T elapsed: ' num2str(tElapsed) ' seconds.'])
end
runInfo.init=init;runInfo.mu=mu;runInfo.sigma=sigma;
runInfo.sigObs=sigObs;runInfo.nsamples=nsamples;runInfo.burn=burn;runInfo.thin=thin;
runInfo.glm=glm;runInfo.iF=iF;
if nargin > 0
    save(saveStr,'chainAll','y_hat_full','runInfo')
end
end


