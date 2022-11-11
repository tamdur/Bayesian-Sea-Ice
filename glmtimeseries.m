function y = glmtimeseries(glm,theta,t)
%GLMTIMESERIES Return glm timeseries given a set of parameter
%specifications
%   INPUTS:
%   glm:    function specification of logistic form to be used
%   theta:  NxK array of parameter estimates, for N draws and k parameters
%   t:      Vector Mx1 of times at which function is evaluated
%
%   OUTPUTS:
%   y:      NxM array of evaluated function.
%
%   Ted Amdur
%   11/10/22

N=size(theta,1);
y=zeros(N,length(t));
for ii=1:size(theta,1)
    y(ii,:) = glm(theta(ii,:),t);
end

end

