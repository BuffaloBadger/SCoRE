function [k0, k0_ci, E, E_ci, r_squared] = Arrhenius_parameters(k,T,R)
%Arrhenius_parameters Fits Arrhenius expression to k vs. T data
%   T must be in K or °R and R must use the same temperature units
%   the energy and moles units for E will be the same as those in R

    % get the number of data
    n_data = length(k);

    % define x and y in the linearized Arrhenius expression
    x = -1./(R*T);
    y = log(k);

    % add a column of ones to x
    x = [x, ones(n_data,1)];

    % fit a linear model to the data 
    [beta, beta_ci, ~, ~, stats] = regress(y,x);

    % extract the Arrhenius parameters
    k0 = exp(beta(2));
    E = beta(1);
    k0_ci = exp(beta_ci(2,:));
    E_ci = beta_ci(1,:);

    % extract the coefficient of determination
    r_squared = stats(1);
end