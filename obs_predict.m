
clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Customize the model run by modifying the following parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define interval examined, slice sampler parameters
plotFlag=1;
cmipBackground=1;
yrsObserved=(1:42)'; %Years in the dataset counted as observations 
nsamples=10000;burn=1000;thin=5;%Number of draws; number of draws to first burn;thinning parameter
iF=1; %inflation factor
sigObs=0.6511.*iF; %std from estimate_obs_std4
meanObs=-0.2681; %mean from estimate_obs_std4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define the model
glm  = @(b,t) b(1).*(1-1./(1+exp(-b(2)*(t-b(3)))).^(1/b(4))); % SIA sigmoid

%Observation values
load sia_obs.mat
y=squeeze(sum(sia_obs,[1 2],'omitnan'))./1e6;
t=yrsObserved;
tAll=(1:122)';

%Draw Priors
load model_full_priorsb_22_10_17.mat


p1=fitdist(p1_priors(:),'Lognormal');
p2=fitdist(p2_priors(:),'Lognormal');
p3=fitdist(p3_priors(:),'Normal');
p4=fitdist(p4_priors(:),'Half Normal');
mu1=p1.mean;mu2=p2.mean;mu3=p3.mean;mu4=p4.mean;
sig1=p1.std;sig2=p2.std;sig3=p3.std;sig4=p4.std;
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
logpost=   @(b) sum(log(normpdf(y-glm(b,t),meanObs,sigObs)))+log(prior1(b(1)))+...
    log(prior2(b(2)))+log(prior3(b(3)))+log(prior4(b(4)));
logprior=  @(b) log(prior1(b(1)))+...
    log(prior2(b(2)))+log(prior3(b(3)))+log(prior4(b(4)));

%Run slice sampler
[chain,neval]=slicesample(init,nsamples,'logpdf',logpost,'thin',thin,'burnin',burn);


for ii = 1:size(chain,1)
    yhat(ii,:) = glm(chain(ii,:),tAll);
end
yint=prctile(yhat,[7,93]);


if plotFlag
    figure
    if cmipBackground
        %Load CMIP6 model data
        load A20.mat
        
        % creating array of model values
        nRuns=79;
        t = yrsObserved;
        tAll=(1:122)';
        ct = 1;
        y_full=zeros(length(t),nRuns);
        y_full=zeros(length(tAll),nRuns);
        for ct_mod = 1:size(A20,2)
            for ct_em = 1:size(A20(ct_mod).X,2)
                y_full(:,ct) = A20(ct_mod).X(tAll,ct_em);
                h(1)=plot(tAll+1978,y_full(:,ct),'LineWidth',0.5,'Color',[0.9 0.9 0.9]);
                hold on
                ct = ct+1;
            end
        end
        hold on
    end
    y50=prctile(yhat,50);
    plot(tAll+1978,yint','Color','k');
    hold on
    h(2)=plot(tAll+1978,y50,'Color','r');
    hold on
    h(3)=scatter(t+1978,y,'MarkerFaceColor','k','MarkerEdgeColor','k');
    xlabel('Year')
    ylabel('SIA')
    set(gca,'FontSize',16)
    ylabel('Sea Ice Area (millions of sq. km)')
    xlim([1978 2100])
    legend(h,'Model runs','Prediction','Observations')
    legend boxoff
end



