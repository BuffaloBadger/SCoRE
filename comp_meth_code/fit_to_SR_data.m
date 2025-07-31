function [beta, betaCI, rSquared] = fit_to_SR_data(betaGuess, x...
    , yMeas, predY, useRelErr)
%fit_to_SR_data fit a model to single response data
%   7/10/24
%
%   betaGuess = column vector containing guesses for the parameter values
%   x = set variable matrix; x(i,j) is the value of set variable j in
%       experiment i
%   yMeas = measured response vector
%   predY = function handle for a function that calculates the model-
%       predicted response vector given the parameter vector and the set
%       variable matrix
%   useRelErr = flag to use relative errors if true
%   beta = column vector of parameter values
%   betaCI = matrix with one row per parameter with the lower limit of 
%       the 95% confidence interval in the first column and the upper 
%       limit in the second column
%   rSquared = coefficient of determination for the fitted model
%
    
    % Calculate the maximum likelihood parameter estimates and their 95%
    % confidence intervals
    if useRelErr
        weight = 1./yMeas;
    else
        weight = ones(size(yMeas,1),1);
    end
    [beta,~,~,sigma,~] = nlinfit(x,yMeas,predY,betaGuess,'Weights',weight);
    eps = (yMeas - predY(beta,x));
    betaCI = nlparci(beta,eps,'covar',sigma);
    
    % Calculate the coefficient of determination
    rSquared = 1.0 - (sum(eps.^2))/(sum((yMeas-mean(yMeas)).^2));
    
end % of fitSR

