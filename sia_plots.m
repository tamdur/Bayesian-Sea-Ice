%Make plots for SIA manuscript. Modified from sigmoidslice_obs plotting
%code
% Ted Amdur and Charlotte Dyvik Henke
% 11/10/22



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings for initial plot showing observations and multi-model
% uncertainty
Fig1=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings for plot using CMIP6 output
Fig2=0; %Set to 1 to use cmip6 output, 0 otherwise
modelrun=24; % Set to the index of the model run used
if Fig2
    load model_predictions_22_11_10.mat %Select set of posterior predictions to run
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings for plot using actual observations
Fig3=0; %Set to 1 to use actual remote sensing observations
newPred=0; %Set to 1 to use new prediction, 0 to load past prediction
if newPred
    %Set parameters for slicesampler
    nsamples=10000;burn=1000;thin=5;%Number of draws; number of draws to first burn;thinning parameter
    saveStr='obs_pred_22_11_11.mat';
    [chain,~,yObs]=obs_predict(nsamples,burn,thin,saveStr);
else
    load obs_pred_22_11_11.mat %Insert name of saved prediction one desires to use
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Prepare data to be plotted
t = yrObs;
tAll = (1:122)';
yrObs=(1:42)'; %Years of observations to be used
if Fig2
    nRuns=79; %Number of CMIP6 model runs to be used
    glm=runInfo.glm; %The functional form used by the MCMC model
    chain=squeeze(chainAll(:,:,modelrun));
    ypred = glmtimeseries(glm,chain,tAll);
    
    %Load CMIP6 model observations
    load A20.mat
    % creating array of model values
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
    %Get full fit
    loglike=   @(b) sum(-1.*log(normpdf(y_full(:,modelrun)-glm(b,tAll),0,runInfo.sigObs)));
    [bfit,~] = fminsearch(loglike, mean(chain,1)');
    y_fit = glm(bfit,tAll);
elseif Fig3
    ypred=glmtimeseries(glm,chain,tAll);
end
siaprctpred = prctile(ypred,[7 50 93])'; %Calculate percentiles

% FIGURE
figure2('Position',[10 10 800 1000])
subplot(2,4,(5:8))
%Shade the 95% CI prediction of logistic
shade(tAll(42:end)+1978,siaprctpred(42:end,1),tAll(42:end)+1978,siaprctpred(42:end,3),'color',rgb('LightPink'),'FillType',[2 1]);
hold on
shade(t+1978,siaprctpred(yrObs,1),t+1978,siaprctpred(yrObs,3),'color',rgb('Gray'),'FillType',[2 1]);
hold on
if Fig2
    h(1)=plot(tAll+1978,y_full(:,modelrun),'k'); %the actual simulated values
else
    h(1)=plot(t+1978,y,'k');
end
hold on
h(2)=plot(tAll+1978,siaprctpred(:,2),'color',rgb('HotPink')); %the 50th percentile of Bayesian fit
hold on
if Fig2
    h(3)=plot(tAll+1978,y_fit,'color',rgb('RoyalBlue')); %full-time bayesian fit
end
yline(1,'--k');
axis tight;
xlabel('Year');
ylabel('SIA [10^6 km^2]');
title('(e)');
ax = gca;
ax.TitleHorizontalAlignment = 'left';
if Fig2
legend(h,'observations','forecast median','122-year fit');
else
    legend('observations','training interval','forecast median','forecast 80% ci')
end


% plots histograms of priors + posteriors for each of the parameters
subplot(2,4,1);hold on;
% setting number of bins for histogram so less dense and colour actually
% shows up
nbins_pr=50;
nbins_po=20;
h1 = histogram(p1_priors,'FaceColor',rgb('Grey'));%priors
h1.NumBins= nbins_pr;
h1.Normalization = 'pdf';
h2 = histogram(chain(:,1),'FaceColor',rgb('HotPink')); %posterior
h2.NumBins = nbins_po;
h2.Normalization = 'pdf';
if Fig2
    xline(bfit(1),'LineWidth',1.5,'color',rgb('RoyalBlue'));
end
title('(a)');
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ylabel('distribution');
xlabel('A_{0} [10^6 km^2]');
xlim([0 20])
legend('prior','posterior','full-fit');

subplot(2,4,2); hold on;
h1 = histogram(p2_priors,'FaceColor',rgb('Grey'));%priors
h1.NumBins= nbins_pr;
h1.Normalization = 'pdf';
h2 = histogram(chain(:,2),'FaceColor',rgb('HotPink')); %posterior
h2.NumBins = nbins_po;
h2.Normalization = 'pdf';
if Fig2
    xline(bfit(2),'LineWidth',1.5,'color',rgb('RoyalBlue'));
end
title('(b)');
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xlim([0 0.5])
ylabel('distribution');
xlabel('alpha [units]');


subplot(2,4,3); hold on;
h1 = histogram(p3_priors+1978,'FaceColor',rgb('Grey'));%priors
h1.NumBins= nbins_pr;
h1.Normalization = 'pdf';
h2 = histogram(chain(:,3)+1978,'FaceColor',rgb('HotPink')); %posterior
h2.NumBins = nbins_po;
h2.Normalization = 'pdf';
if Fig2
    xline(bfit(3)+1978,'LineWidth',1.5,'color',rgb('RoyalBlue'));
end
title('(c)');
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xlim([1950 2060])
ylabel('distribution');
xlabel('t_{0} [year]');



subplot(2,4,4); hold on;
h1 = histogram(p4_priors,'FaceColor',rgb('Grey'));%priors
h1.NumBins= nbins_pr;
h1.Normalization = 'pdf';
h2 = histogram(chain(:,4),'FaceColor',rgb('HotPink')); %posterior
h2.NumBins = nbins_po;
h2.Normalization = 'pdf';
if Fig2
    xline(bfit(4),'LineWidth',1.5,'color',rgb('RoyalBlue'));
end
title('(d)');
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xlim([0 10])
ylabel('distribution');
xlabel('\nu [unitless]');
    

% saving the figure
FileName = 'SIAobs_bayes_22_11_10';
saveas(gcf,['plots/' FileName],'png');
