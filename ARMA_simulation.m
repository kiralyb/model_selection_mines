function ARMA_simulation
% ARMA_SIMULATION     Using Bayesian information criterion (BIC) and Akaike
% information criterion (AIC) to choose the best fitting
% autoregressive-moving-average (ARMA) model. 

% Local field potential data is simulated with an AMRA (p=2,q=1) process
% and ARMA models are fitted on the data with different p and q values
% in the range 1 to max_order. BIC is calculated to find the best matching
% model which is then used to predict the data.

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   18-May-2022

% Define parameters
max_order = 4; % maximum p and q values used for model fitting
time_length = 80; % length of the time series data

% Simulate ARMA(2,1) data
Model = arima('Constant',1,'AR',{0.5 -0.20},'MA',0.4,'Variance',0.1);
rng(15)
y = simulate(Model,time_length);

% Fit ARMA models on the data with all combinations of p and q =
% 1,...,max_order and calculate the corresponding loglikelyhood values
LgLikelyHood = zeros(max_order);
PQ = zeros(max_order);
for p = 1:max_order
    for q = 1:max_order
        Mdl = arima(p,0,q);
        [~,~,logL] = estimate(Mdl,y,'Display','off');
        LgLikelyHood(p,q) = logL;
        PQ(p,q) = p + q;
    end
end

% Calculate and plot Bayesian information criteria (BIC) for each fitted
% model
LgLikelyHood = reshape(LgLikelyHood,max_order^2,1);
PQ = reshape(PQ,max_order^2,1);
[aic,bic] = aicbic(LgLikelyHood,PQ + 1,time_length);
figure
imagesc(reshape(bic,max_order,max_order))
title('BIC')
colorbar
figure
imagesc(reshape(aic,max_order,max_order))
title('AIC')
colorbar

% Find the p and 1 order of the best fitting model (minimal BIC)
[~,p] = min(min(reshape(bic,max_order,max_order),[],2));
[~,q] = min(min(reshape(bic,max_order,max_order),[],1));


% fit the best fitting ARMA (p,q) model and use it for prediction
Mdl = arima(p,0,q);
[EstMdl,~,~] = estimate(Mdl,y,'Display','off');
prediction = zeros(1,time_length);
for i = p + q + 1:time_length
    prediction(i) = forecast(EstMdl,1,'Y0',y(i - p - q:i - 1));
end

% plot the simulated data and the prediction of the best matching model
figure
plot(1:length(y),y)
hold on
plot(p + q + 1:length(y),prediction(p + q + 1:end))
