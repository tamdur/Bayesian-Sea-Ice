function [p1_priors,p2_priors,p3_priors,p4_priors,runInfo]= sigmoid_params(A20,saveStr)

    % creating a vector for models
    t = (1:122)';
    ct = 1;
    for ct_mod = 1:20
        for ct_em = 1:size(A20(ct_mod).X,2)
            y(:,ct) = A20(ct_mod).X(:,ct_em);
            ct = ct+1;
        end
    end
    
    
    %Sampler parameters
    sigObs = 0.6511; %From estimate_obs_std4.m

    %Number of draws; number of draws to first burn;thinning parameter 
    nsamples=10000;burn=1000;thin=5;    

    % rough priors
    yhat=zeros(length(t),nsamples);
    p1_priors=zeros(nsamples,size(y,2));
    p2_priors=zeros(nsamples,size(y,2));
    p3_priors=zeros(nsamples,size(y,2));
    p4_priors=zeros(nsamples,size(y,2));
    
    tElapsed=0;
    for ii = 1:79
        tic
    % ====================================================================== %
      % set initial values, mean, sigma.  

                  %b1 - beta    %b2 - alpha    %b3 - t
        init    = [6            0.5          70       10];
                  %mu1          %mu2           %mu3
        mu      = [6            0          70       0];
                  %sigma1       %sigma2        %sigma3  
        sigma   = [10           0.5          70       3];

        %Then, define the prior pdfs for the parameters
        disttypes=["HalfNormal";"Normal";"HalfNormal";"Normal"];
        halfnormpdf2=makedist('HalfNormal','mu',mu(2),'sigma',sigma(2));
        halfnormpdf4=makedist('HalfNormal','mu',mu(4),'sigma',sigma(4));
        prior1= @(b1) normpdf(b1,mu(1),sigma(1));
        prior2= @(b2) pdf(halfnormpdf2,b2);
        prior3= @(b3) normpdf(b3,mu(3),sigma(3));
        prior4= @(b4) pdf(halfnormpdf4,b4);

        % Define the model (with 3 free parameters)
        glm  = @(b,t) b(1).*(1-1./(1+exp(-b(2)*(t-b(3)))).^(1/b(4))); % SIA sigmoid


        %Write function for likelihood
        loglik=    @(b) sum(log(normpdf(y(:,ii)-glm(b,t),0,sigObs)))+log(prior2(b(2)))+log(prior4(b(4)));


        %Run slice sampler
        clear yhat
        tic
        [chain_p,~]=slicesample(init,nsamples,'logpdf',loglik,'thin',thin,'burnin',burn);
        for iR = 1:size(chain_p,1)
            yhat(:,iR) = glm(chain_p(iR,:),t);
        end
        toc %Display elapsed time spent sampling

        p1_priors(:,ii) = chain_p(:,1);
        p2_priors(:,ii) = chain_p(:,2);
        p3_priors(:,ii) = chain_p(:,3);
        p4_priors(:,ii) = chain_p(:,4);
        
        nElapsed=toc;
        tElapsed=tElapsed+nElapsed;
        disp(['Run ' num2str(ii) ' posterior complete. T elapsed: ' num2str(tElapsed) ' seconds.'])
    end
    runInfo.init=init;runInfo.mu=mu;runInfo.sigma=sigma;
    runInfo.sigObs=sigObs;runInfo.nsamples=nsamples;runInfo.burn=burn;runInfo.thin=thin;
    runInfo.glm=glm;runInfo.disttypes=disttypes;
    if nargin>1
        save(saveStr,'p1_priors','p2_priors','p3_priors','p4_priors','runInfo')
    end
end



