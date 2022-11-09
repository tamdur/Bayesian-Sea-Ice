%Estimate the standard deviation of the observations from the long-term
%sigmoidal trend by examining average residuals from CMIP6 model runs

clearvars
% TESTING THE BAYESIAN SAMPLER ON EACH OF THE MODELS

%Perform a leave-one-out prediction of

% Define the model
glm  = @(b,t) b(1).*((1-1./(1+exp(-b(2)*(t-b(3)))))); % SIA sigmoid

%Get the maximum likelihood fit for each model
load model_predictions_22_10_18.mat
load A20.mat
ct = 1;
for ct_mod = 1:size(A20,2)
    for ct_em = 1:size(A20(ct_mod).X,2)
        y(:,ct) = A20(ct_mod).X(1:42,ct_em);
        y2(:,ct) = A20(ct_mod).X(:,ct_em);
        ct = ct+1;
    end
end
t = (1:42)';
t2 = (1:122)';
% % storing most likely 1st round priors
% for ct_79 = 1:79
%     pd1=fitdist(p1_priors(:,ct_79),'Normal');
%     mlp1_full(ct_79)=pd1.mean;
%     pd2=fitdist(p2_priors(:,ct_79),'Normal');
%     mlp2_full(ct_79)=pd2.mean;
%     pd3=fitdist(p3_priors(:,ct_79),'Normal');
%     mlp3_full(ct_79)=pd3.mean;
% 
%     %Generate maximum likelihood logistic fit to each realization
%     y_hat_full(:,ct_79) = glm([mlp1_full(ct_79) mlp2_full(ct_79) mlp3_full(ct_79)],t2);
% end

%Get a distribution of residuals
residuals=y2(1:42,:)-y_hat_full(1:42,:);
%standard distance from zero:
res=sqrt(sum(residuals(:).^2)./length(residuals(:))); %yields 0.6511
%std(residuals) yields 0.4826



