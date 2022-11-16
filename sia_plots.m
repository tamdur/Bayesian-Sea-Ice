%Make plots for SIA manuscript. Modified from sigmoidslice_obs plotting
%code
% Ted Amdur and Charlotte Dyvik Henke
% 11/10/22

clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings for initial plot showing observations and multi-model
% uncertainty
Fig1=0;
if Fig1
saveStr='SIAFig1_22_11_15';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings for plotting a single CMIP6 output run
Fig2=0; %Set to 1 to use cmip6 output, 0 otherwise
modelrun=19; % Set to the index of the model run used
if Fig2
saveStr='SIAFig2_22_11_15';
load model_predictions_22_11_14b.mat %Select set of posterior predictions to run
load model_full_priors_22_11_14.mat
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings for plotting a panel of all the CMIP6 runs
Fig2b=0;
if Fig2b
    saveStr='SIAFig2b_22_11_15';
    load model_predictions_22_11_14b.mat %Select set of posterior predictions to run
    load model_full_priors_22_11_14.mat
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings for plot using actual observations
Fig3=0; %Set to 1 to use actual remote sensing observations
newPred=0; %Set to 1 to use new prediction, 0 to load past prediction
if Fig3
saveStr='SIAFig3_22_11_16';
load model_full_priors_22_11_14.mat %Priors to plot
end
if newPred
    %Set parameters for slicesampler
    nsamples=10000;burn=1000;thin=5;%Number of draws; number of draws to first burn;thinning parameter
    saveStr='obs_pred_22_11_15.mat';
    [chain,~,yObs]=obs_predict(nsamples,burn,thin,saveStr);
else
    load obs_pred_22_11_16long.mat %Insert name of saved prediction one desires to use
    glm=obsInfo.glm;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings for plot of sensitivity to parameters
Fig4=1; %Set to 1 to use actual remote sensing observations
newPred=0; %Set to 1 to use new prediction, 0 to load past prediction
fSize=16; %Font Size
if Fig4
saveStr='SIAFig4_22_11_16';
load model_full_priors_22_11_14.mat %Priors to plot
end
if newPred
    %Set parameters for slicesampler
    nsamples=10000;burn=1000;thin=5;%Number of draws; number of draws to first burn;thinning parameter
    saveStr='obs_pred_22_11_15.mat';
    [chain,~,yObs]=obs_predict(nsamples,burn,thin,saveStr);
else
    load obs_pred_22_11_16long.mat %Insert name of saved prediction one desires to use
    glm=obsInfo.glm;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Prepare data to be plotted
yrObs=(1:42)'; %Years of observations to be used
t = yrObs;
tAll = (1:122)';
stYr=1978;

%Script for Fig. 1
if Fig1
    load sia_obs.mat;
    load A20.mat;
    
    y_obs = squeeze(sum(sia_obs, [1 2],'omitnan')./1e6);
    ct = 1;
    for ct_mod = 1:20
        for ct_em = 1:size(A20(ct_mod).X,2)
            y(:,ct) = A20(ct_mod).X(1:42,ct_em);
            y2(:,ct) = A20(ct_mod).X(:,ct_em);
            ct = ct+1;
        end
    end
    % ======================================================================
    % Get priors for each of the 79 realisations
    load model_full_priorsb_22_11_10.mat
    
    
    % storing 1st round priors, finding most frequent value
    glm=runInfo.glm; %The functional form used by the MCMC model
    for ct_79 = 1:79
        %Get full fit
        mle=[mean(p1_priors(:,ct_79)) mean(p2_priors(:,ct_79)) mean(p3_priors(:,ct_79)) mean(p4_priors(:,ct_79))];
        y_full(:,ct_79) = glm(mle,tAll);
    end
    siamodsprct=prctile(y_full',[2.5 50 97.5]);
    
    % creating index for first ensemble member for each model
    for ct_mod = 1:20
        em1(ct_mod) = size(A20(ct_mod).X,2);
    end
    em1= cumsum(em1);
    
    % plotting figure 1
    figure; hold on;
    plot(t+stYr,y_obs,'color',rgb('CornflowerBlue'),'LineWidth',1.3);%observations
    plot(tAll+stYr,siamodsprct(2,:),'color',rgb('Crimson'),'LineWidth',1.3);%multimodel median
    shade(tAll+stYr,siamodsprct(1,:),tAll+stYr,siamodsprct(3,:),'color',rgb('LightPink'),'FillType',[2 1]);
    for ct =em1
        plot(tAll+stYr,y2(:,ct),'color',[rgb('Black'),0.15]); %actual realisations
        plot(tAll+stYr,y_full(:,ct),':','color',rgb('Black')); %logistic fits
        yline(1,'--k');%threshold
        %     plot(t2+stYr,siamodsprct(:,2),'color',rgb('Red'),'LineWidth',1.5);%model median again
        axis tight;
        xlabel('Year');
        ylabel('Sea ice area [10^6 km^2]');
        yline(1,'--k');
    end
    plot(t+stYr,y_obs,'color',rgb('CornflowerBlue'),'LineWidth',1.3);%observations again
    plot(tAll+stYr,siamodsprct(2,:),'color',rgb('Crimson'),'LineWidth',1.3);%multimodel median again
    legend('observations','multi-model median','','','multi-model 95% c.i.','simulations','logistic fits','sea-ice-free threshold');
    
    % saving the figures
    % saving the figure
    saveas(gcf,['plots/' saveStr],'png');
end
if Fig2b
    glm=obsInfo.glm;
    figure
    for mn=1:79
        full_fit_theta=[median(p1_priors(:,mn)) median(p2_priors(:,mn)) median(p3_priors(:,mn)) median(p4_priors(:,mn))]; 
        yhat = glmtimeseries(glm,squeeze(chainAll(:,:,mn)),tAll);
        yfull=glmtimeseries(glm,full_fit_theta,tAll);
        yint=prctile(yhat,[5 95]);
        subplot(10,8,mn)
        plot(stYr+tAll,yint','Color','k');
        hold on; 
        plot(stYr+tAll,yfull,'Color','r');
        hold on;
        title(['Run ' num2str(mn)])
    end 
    saveas(gcf,['plots/' saveStr],'png');
end
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
    bfit = mean(chain,1)';
    y_fit = y_hat_full(:,modelrun);
    siaprctpred = prctile(ypred,[7 50 93])'; %Calculate percentiles
    % FIGURE
    figure2('Position',[10 10 800 1000])
    subplot(2,4,(5:8))
    %Shade the 95% CI prediction of logistic
    shade(tAll(42:end)+stYr,siaprctpred(42:end,1),tAll(42:end)+stYr,siaprctpred(42:end,3),'color',rgb('LightPink'),'FillType',[2 1]);
    hold on
    shade(t+stYr,siaprctpred(yrObs,1),t+stYr,siaprctpred(yrObs,3),'color',rgb('Gray'),'FillType',[2 1]);
    hold on
    h(1)=plot(tAll+stYr,y_full(:,modelrun),'k'); %the actual simulated values
    hold on
    h(2)=plot(tAll+stYr,siaprctpred(:,2),'color',rgb('HotPink')); %the 50th percentile of Bayesian fit
    hold on
    h(3)=plot(tAll+stYr,y_fit,'color',rgb('RoyalBlue')); %full-time bayesian fit
    yline(1,'--k');
    axis tight;
    xlabel('Year');
    ylabel('SIA [10^6 km^2]');
    title('(e)');
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    legend(h,'observations','forecast median','122-year fit');
    
    
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
    xline(bfit(1),'LineWidth',1.5,'color',rgb('RoyalBlue'));
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
    xline(bfit(2),'LineWidth',1.5,'color',rgb('RoyalBlue'));
    title('(b)');
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    xlim([0 0.5])
    ylabel('distribution');
    xlabel('alpha [units]');
    
    
    subplot(2,4,3); hold on;
    h1 = histogram(p3_priors+stYr,'FaceColor',rgb('Grey'));%priors
    h1.NumBins= nbins_pr;
    h1.Normalization = 'pdf';
    h2 = histogram(chain(:,3)+stYr,'FaceColor',rgb('HotPink')); %posterior
    h2.NumBins = nbins_po;
    h2.Normalization = 'pdf';
    xline(bfit(3)+stYr,'LineWidth',1.5,'color',rgb('RoyalBlue'));
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
    xline(bfit(4),'LineWidth',1.5,'color',rgb('RoyalBlue'));
    title('(d)');
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    xlim([0 10])
    ylabel('distribution');
    xlabel('\nu [unitless]');

    % saving the figure
    saveas(gcf,['plots/' saveStr],'png');

end
if Fig3
        ypred=glmtimeseries(glm,chain,tAll);
    siaprctpred = prctile(ypred,[7 50 93])'; %Calculate percentiles
    
    % FIGURE
    figure2('Position',[10 10 800 1000])
    subplot(2,4,(5:8))
    %Shade the 95% CI prediction of logistic
    shade(tAll(42:end)+stYr,siaprctpred(42:end,1),tAll(42:end)+stYr,siaprctpred(42:end,3),'color',rgb('LightPink'),'FillType',[2 1]);
    hold on
    shade(t+stYr,siaprctpred(yrObs,1),t+stYr,siaprctpred(yrObs,3),'color',rgb('Gray'),'FillType',[2 1]);
    hold on

        h(1)=plot(t+stYr,yObs,'k');
    hold on
    h(2)=plot(tAll+stYr,siaprctpred(:,2),'color',rgb('HotPink')); %the 50th percentile of Bayesian fit
    hold on
    yline(1,'--k');
    axis tight;
    xlabel('Year');
    ylabel('SIA [10^6 km^2]');
    title('(e)');
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
        legend('observations','training interval','forecast median','forecast 80% ci')
    
    
    % plots histograms of priors + posteriors for each of the parameters
    subplot(2,4,1);hold on;
    % setting number of bins for histogram so less dense and colour actually
    % shows up
    nbins_pr=50;
    nbins_po=50;
    h1 = histogram(p1_priors,'FaceColor',rgb('Grey'));%priors
    h1.NumBins= nbins_pr;
    h1.Normalization = 'pdf';
    h2 = histogram(chain(:,1),'FaceColor',rgb('HotPink')); %posterior
    h2.NumBins = nbins_po;
    h2.Normalization = 'pdf';
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
    title('(b)');
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    %set(ax,'xscale','log')
    xlim([0 2])
    ylabel('distribution');
    xlabel('alpha [units]');
    
    
    subplot(2,4,3); hold on;
    h1 = histogram(p3_priors+stYr,'FaceColor',rgb('Grey'));%priors
    h1.NumBins= nbins_pr;
    h1.Normalization = 'pdf';
    h2 = histogram(chain(:,3)+stYr,'FaceColor',rgb('HotPink')); %posterior
    h2.NumBins = nbins_po;
    h2.Normalization = 'pdf';
    title('(c)');
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    xlim([1970 2080])
    ylabel('distribution');
    xlabel('t_{0} [year]');
    
    
    
    subplot(2,4,4); hold on;
    h1 = histogram(p4_priors,'FaceColor',rgb('Grey'));%priors
    h1.NumBins= nbins_pr;
    h1.Normalization = 'pdf';
    h2 = histogram(chain(:,4),'FaceColor',rgb('HotPink')); %posterior
    h2.NumBins = nbins_po;
    h2.Normalization = 'pdf';
    title('(d)');
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    %set(ax,'xscale','log')
    xlim([0 50])
    ylabel('distribution');
    xlabel('\nu [unitless]');
    
    
    % saving the figure
    saveas(gcf,['plots/' saveStr],'png');
end
if Fig4
    yhat=glmtimeseries(glm,chain,tAll);
    %Find closest year for each draw
    [~,yr]=min(abs(yhat-1),[],2);
    priorsAll=[p1_priors(:) p2_priors(:) p3_priors(:) p4_priors(:)];
    yhat_all=glmtimeseries(glm,priorsAll(1:10:end,:),tAll);
    [~,yr_all]=min(abs(yhat_all-1),[],2);
    
    ind=1;
    for ii=min(yr):max(yr)
        condI(ind).iYr=find(yr==ii);
        condI(ind).year=ii;
        condI(ind).conRange1=prctile(chain(condI(ind).iYr,1),[10 50 90]);
        condI(ind).conRange2=prctile(chain(condI(ind).iYr,2),[10 50 90]);
        condI(ind).conRange3=prctile(chain(condI(ind).iYr,3),[10 50 90]);
        condI(ind).conRange4=prctile(chain(condI(ind).iYr,4),[10 50 90]);
        ind=ind+1;
    end
    
    %Start by plotting just the first parameter, the max asymptote 
    xEdges=linspace(min(yr)-1,max(yr)+1,max(yr)-min(yr)+3); %Add an edge before an after so line varies from zero
    xPts=diff(xEdges)+xEdges(1:end-1);
    x=histcounts(yr,xEdges);
    xprior=histcounts(yr_all,xEdges);
    figure2('Position',[10 10 1000 1000])
    subplot(2,2,1)
    pSel=1; %Parameter being plotted
    yEdges=linspace(min(chain(:,pSel)),max(chain(:,pSel)),30);
    yEdges=[yEdges(1)-range(yEdges)./30 yEdges yEdges(end)+range(yEdges)./30];
    yPts=diff(yEdges)+yEdges(1:end-1);
    y=histcounts(chain(:,pSel),yEdges);
    yprior=histcounts(p1_priors(:),yEdges);
    
    for ii=1:length(condI)
        xYr(ii)=condI(ii).year+stYr;
        y10(ii)=condI(ii).conRange1(1);
        y50(ii)=condI(ii).conRange1(2);
        y90(ii)=condI(ii).conRange1(3);
    end
    plot(xYr,y10,'r--')
    hold on
    plot(xYr,y50,'Color','k','LineWidth',2)
    hold on
    plot(xYr,y90,'r--')
    hold on
    xy=x*(range(chain(:,pSel))./(10*max(x)))+min(yPts)-range(chain(:,pSel))./10;
    xyprior=xprior*(range(chain(:,pSel))./(10*max(xprior)))+min(yPts)-range(chain(:,pSel))./10;
    plot(xPts+stYr,xyprior,'b--')
    hold on
    plot(xPts+stYr,xy,'Color','b')
    hold on
    plot(yprior*(range(yr)./(0.1*max(yprior)*max(yr)))+min(yr)+stYr-5,yPts,'b--');
    hold on
    plot(y*(range(yr)./(0.1*max(y)*max(yr)))+min(yr)+stYr-5,yPts,'Color','b');
    xlim([2019 2060])
    ylim([min(xy) max(yPts)])
    ylabel('A_{0} [10^6 km^2]');
    xlabel('Year of Ice-Free Arctic')
    set(gca,'FontSize',fSize)
    
    subplot(2,2,2)
    pSel=2; %Parameter being plotted
    yEdges=linspace(min(chain(:,pSel)),max(chain(:,pSel)),30);
    yEdges=[yEdges(1)-range(yEdges)./30 yEdges yEdges(end)+range(yEdges)./30];
    yPts=diff(yEdges)+yEdges(1:end-1);
    y=histcounts(chain(:,pSel),yEdges);
    yprior=histcounts(p2_priors(:),yEdges);
    
    for ii=1:length(condI)
        xYr(ii)=condI(ii).year+stYr;
        y10(ii)=condI(ii).conRange2(1);
        y50(ii)=condI(ii).conRange2(2);
        y90(ii)=condI(ii).conRange2(3);
    end
    plot(xYr,y10,'r--')
    hold on
    plot(xYr,y50,'Color','k','LineWidth',2)
    hold on
    plot(xYr,y90,'r--')
    hold on
    xy=x*(range(chain(:,pSel))./(10*max(x)))+min(yPts)-range(chain(:,pSel))./10;
    xyprior=xprior*(range(chain(:,pSel))./(10*max(xprior)))+min(yPts)-range(chain(:,pSel))./10;
    plot(xPts+stYr,xyprior,'b--')
    hold on
    plot(xPts+stYr,xy,'Color','b')
    hold on
    plot(yprior*(range(yr)./(0.1*max(yprior)*max(yr)))+min(yr)+stYr-5,yPts,'b--');
    hold on
    plot(y*(range(yr)./(0.1*max(y)*max(yr)))+min(yr)+stYr-5,yPts,'Color','b');
    xlim([2019 2060])
    ylim([min(xy) max(yPts)])
    xlabel('Year of Ice-Free Arctic')
    ylabel('\alpha');
    set(gca,'FontSize',fSize)
    
    subplot(2,2,3)
    pSel=3; %Parameter being plotted
    yEdges=linspace(min(chain(:,pSel)),max(chain(:,pSel)),30);
    yEdges=[yEdges(1)-range(yEdges)./30 yEdges yEdges(end)+range(yEdges)./30];
    yPts=diff(yEdges)+yEdges(1:end-1)+stYr;
    y=histcounts(chain(:,pSel)+stYr,yEdges+stYr);
    yprior=histcounts(p3_priors(:),yEdges);
    
    for ii=1:length(condI)
        xYr(ii)=condI(ii).year+stYr;
        y10(ii)=condI(ii).conRange3(1)+stYr;
        y50(ii)=condI(ii).conRange3(2)+stYr;
        y90(ii)=condI(ii).conRange3(3)+stYr;
    end
    plot(xYr,y10,'r--')
    hold on
    plot(xYr,y50,'Color','k','LineWidth',2)
    hold on
    plot(xYr,y90,'r--')
    hold on
    xy=x*(range(chain(:,pSel))./(10*max(x)))+min(yPts);
    xyprior=xprior*(range(chain(:,pSel))./(10*max(xprior)))+min(yPts);
    plot(xPts+stYr,xyprior,'b--')
    hold on
    plot(xPts+stYr,xy,'Color','b')
    hold on
    plot(yprior*(range(yr)./(0.1*max(yprior)*max(yr)))+min(yr)+stYr-5,yPts,'b--');
    hold on
    plot(y*(range(yr)./(0.1*max(y)*max(yr)))+min(yr)+stYr-5,yPts,'Color','b');
    xlim([2019 2060])
    ylim([min(xy) max(yPts)])
    xlabel('Year of Ice-Free Arctic')
    ylabel('t_{0} [year]');
    set(gca,'FontSize',fSize)
    
    subplot(2,2,4)
    pSel=4; %Parameter being plotted
    yEdges=linspace(min(chain(:,pSel)),max(chain(:,pSel)),30);
    yEdges=[yEdges(1)-range(yEdges)./30 yEdges yEdges(end)+range(yEdges)./30];
    yPts=diff(yEdges)+yEdges(1:end-1);
    y=histcounts(chain(:,pSel),yEdges);
    yprior=histcounts(p4_priors(:),yEdges);
    
    for ii=1:length(condI)
        xYr(ii)=condI(ii).year+stYr;
        y10(ii)=condI(ii).conRange4(1);
        y50(ii)=condI(ii).conRange4(2);
        y90(ii)=condI(ii).conRange4(3);
    end
    plot(xYr,y10,'r--')
    hold on
    plot(xYr,y50,'Color','k','LineWidth',2)
    hold on
    plot(xYr,y90,'r--')
    hold on
    xy=x*(range(chain(:,pSel))./(10*max(x)))+min(yPts)-range(chain(:,pSel))./10;
    xyprior=xprior*(range(chain(:,pSel))./(10*max(xprior)))+min(yPts)-range(chain(:,pSel))./10;
    plot(xPts+stYr,xyprior,'b--')
    hold on
    plot(xPts+stYr,xy,'Color','b')
    hold on
    plot(yprior*(range(yr)./(0.1*max(yprior)*max(yr)))+min(yr)+stYr-5,yPts,'b--');
    hold on
    plot(y*(range(yr)./(0.1*max(y)*max(yr)))+min(yr)+stYr-5,yPts,'Color','b');
    xlim([2019 2060])
    ylim([min(xy) max(yPts)])
    xlabel('Year of Ice-Free Arctic')
    ylabel('\nu')
    set(gca,'FontSize',fSize)
    
    % saving the figure
    saveas(gcf,['plots/' saveStr],'png');
    
end
