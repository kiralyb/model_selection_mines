function tuning_model_sim (gain_factor,additive_factor,simnum)
%TUNING_MODEL_SIM Model fitting of simulated tuning curve change data.
%   TUNING_MODEL_SIM (GAIN_FACTOR,ADDITIVE_FACTOR,SIMNUM) % Simulating
%   neuronal response (tuning curve) change data with a gain model y =
%   ax(s) + b + n (mixture of multiplicative and additive models) where
%   y = x(s) is the baseline tuning curve (bell curve), a is the gain
%   factor,b is additive gain component, n is a Gaussian noise term. An
%   additive and a multiplicative model is fitted on the simulated data and
%   Akaike information criterion (AIC) is calculated to find the better
%   fitting model. The simulation is repeated n times and the mean AIC
%   values are calculated along with the AIC difference distribution. 

%   Required input arguments:
%       GAIN_FACTOR: multiplicative factor applied for data simulation
%       ADDITIVE_FACTOR: additive factor applied for data simulation 
%       SIMNUM: number of simulations

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   18-May-2022


% define parameters
wind = 1:13; % x coordiantes
mu = 7.8; % mean of the baseline Gaussian
sigma = 1.3; % sigma^2 is the variance of the baseline Gaussian
const = 0.5; % additive constant shift // applied to have only positive data
noise_snr = 18; % noise level (additive Gaussian noise) of the simulated data
if simnum > 1
    figvis = 'off';
else
    figvis = 'on';
end

% perform tuning gain model simulation simnum times
Diff = zeros(1,simnum);
AIC = zeros(1,simnum);
for i=1:simnum
    [Diff(i),AIC(i,:)] = tuning_curve_fit(wind,mu,sigma,const,gain_factor,additive_factor,noise_snr,figvis);    
end

if simnum > 1
% plot the histogram of AIC differences
figure
histogram(Diff)

%plot mean AICs with error bars (standard deviation)
figure
bar([1,2],[mean(AIC(:,1)),mean(AIC(:,2))])
hold on
errorbar([1,2],[mean(AIC(:,1)),mean(AIC(:,2))],[std(AIC(:,1)),std(AIC(:,2))])
end



function [diff_,AIC] = tuning_curve_fit(wind,mu,sigma,const,gain_factor,additive_factor,noise_snr,visible)


% simulate baseline tuning data
rng('shuffle')
tuning = normpdf(wind,mu,sigma) + const;

% plot baseline tuning curve
figure('visible',visible)
hold on
fplot(@(x) 1 / sigma / sqrt(2 * pi()) * (exp(-(x - mu)^2/(2 * sigma^2))) + const,[wind(1),wind(end)]);

% simulate and plot multiplicative gain model data
tuning2 = tuning * gain_factor + additive_factor + awgn(tuning,noise_snr);
[xData, yData] = prepareCurveData( [], tuning2 );
scatter(xData, yData,'+','k');

% fit and plot multiplicative gain model tuning curve
f = fittype(['a*(' num2str(1 / sigma / sqrt(2 * pi())) '*(exp(-(x-' ...
    num2str(mu) ')^2/(2*' num2str(sigma) '^2)))+' ...
    num2str(const) ')']);
[fun{1}, gof{1}] = fit( xData, yData, f);
plot(fun{1},'b');

% fit and plot additive gain model tuing curve
f = fittype(['a+(' num2str(1 / sigma / sqrt(2 * pi())) '*(exp(-(x-' ...
    num2str(mu) ')^2/(2*' num2str(sigma) '^2)))+' ...
    num2str(const) ')']);
[fun{2}, gof{2}] = fit( xData, yData, f);
plot(fun{2},'r');

%add legends to the plot
legend('baseline tuning curve','simulated gain data', 'multiplicative gain model', 'additive gain model', 'Location', 'NorthEast' );
YL=get(gca,'YLim');
set(gca,'Ylim',[const, YL(2)]);


% calculate Akaike information criterion (AIC) for the two model
AIC = zeros(1,2);
for i = 1:2
    n = length(tuning);   % number of data points
    k = 1;   % number of model parameters
    R = gof{i}.rmse .^ 2;   % mean square error
    AIC(i) = n * log(R/n) + 2 * k;   % Akaike information criterion
end

% plot AIC and calculate difference
figure('visible',visible)
bar(AIC)
diff_ = AIC(2) - AIC(1);
