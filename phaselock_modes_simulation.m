function phaselock_modes_simulation()
%PHASELOCK_MODES_SIMULATION Simulating phase preference data and fitting
% the mixture of von Mises distributuions

% Phase preference data of neuronal firing referenced to an LFP oscillation
% is simulated as the combination of three wrapped normal distributions.
% The mixture of 1 to max_components von Mises distributions are fitted on the simulated 
% distribution using an expectation maximalization algorithm and AIC, BIC,
% and parametric bootstrap p-value is calculated to find the best fitting
% model. 

%   Required input arguments:
%       MAINPATH: 

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   11-Nov-2021


% parameters 
B_stimnum = 100; % number of bootstrap simulations 
max_components = 4; % maximal number of modes for the fitted distributions  

% simulating the combination of 3 normal distributions 
rng(1)
Phase1 = normrnd(pi,pi/10,60,1);
Phase2 = normrnd(7*pi/4,pi/5,125,1);
Phase3 = normrnd(2.3*pi/4,pi/5,65,1);    % Phase3 = normrnd(pi/2,pi/5,65,1);
Phase = [Phase1;Phase2;Phase3];

% transforming the distribution to a mixtrue of wrapped distributions 
Phase(Phase<0) = Phase(Phase<0)+2*pi;
Phase(Phase>2*pi) = Phase(Phase>2*pi)-2*pi;

% phase histogram plot
figure
rose(Phase)

% fitting mixture of von Mises distributions, and calculating
% AIC,BIC,p-value
[pvalue AIC BIC]=vmcomponents_(Phase,'rad',max_components,B_stimnum)

% plot AIC, BIC, and p in the function of number of modes
figure
plot(AIC)
hold on
plot(BIC)
yyaxis right
plot(pvalue)

function [pvalue AIC BIC] = vmcomponents_(theta,degrad,max_components, B_stimnum)
%VMCOMPONENTS_   Estimates the number of von Mises components.
%   [PVALUE AIC BIC] = VMCOMPONENTS_(THETA,DEGRAD,MAX_COMPONENTS,B) performs statistical tests
%   addressing the question whether the number of modes in the circular
%   distribution of the sample THETA is (i) k or (ii) more than k, from k =
%   1 to MAX_COMPONENTS. The second input argument (DEGRAD) should contain the
%   information whether the sample is measured in degress ('deg') or
%   radians ('rad'). PVALUEs (based on B parametric bootstrap resampling)
%   are returned along with AIC and BIC
%   information criteria for model selection.
%
%   Reference: Fisher NI (1993) Statistical analysis of circular data,
%   Cambridge University Press, Cambridge. pp. 100-102.
%
%   See also WATSONKFIT and WATSONKFIT_EM.

% Input argument check

dbstop if error
rng shuffle
%error(nargchk(3,3,nargin))
if isequal(degrad,'deg')    % convert to radian
    theta = theta / 180 * pi;
end

% Directories
%global DATAPATH
inpdir = 'L:\Balint\NatNeuro_prespective\phaselock'
resdir = 'L:\Balint\NatNeuro_prespective\phaselock'

% Estimate the number of von Mises components
pvalue = zeros(1,max_components);
AIC = zeros(1,max_components);
BIC = zeros(1,max_components);
for k = 1:max_components
    [pvalue(k) AIC(k) BIC(k)] = main(theta,k,inpdir,resdir, B_stimnum);
    save([resdir 'pvalue.mat'],'pvalue','AIC','BIC')
end



% -------------------------------------------------------------------------
function [pvalue AIC BIC] = main(theta,k,inpdir,resdir, B)

% Load parameter estimations (Step 1.)
n = length(theta);

[param,err] = watsonkfit(theta,k,1);
mu=param{1};
kappa=param{2};
p=param{3};

p(end+1) = 1 - sum(p);

F = cell(1,k);    % mixed von Mises cdf
cF = zeros(1,n);
for t = 1:k
    F{t} = vmcdf(theta',mu(t),kappa(t),'rad');
    cF = cF + p(t) * F{t};
end

L = 0;  % log-likelihhod function
for j = 1:k
    L = L + p(j)/(2*pi*besseli(0,kappa(j)))*exp(kappa(j)*cos(theta-mu(j)));
end
Q = sum(log(L));
AIC = 2 * (3 * k - 1) - 2 * Q;
BIC = - 2 * Q + (3 * k - 1) * log(n);
pvalue = 0;

% Goodness-of-fit (Step 2.)
z = sort(cF,'ascend');
T0 = (z - (2 * (1:n) - 1) / (2 * n)) .^ 2;
T = sum(T0) - n * (mean(z) - 0.5) ^ 2 + 1 / (12 * n);   % U^2 (4.35)

% Parametric resample (Step. 3.)

T_star = zeros(1,B);
for bts = 1:B
    disp(['bts = ' num2str(bts)])
    Phat0 = 0;
    Phat = cumsum(p);   % (3.1)
    u = rand(1,n);      % U(0,1) (3.2)
    theta_star = zeros(1,n);
    for t = 1:n         % (3.3)
        j = find([Phat0 Phat]<=u(t),1,'first');
        theta_star(t) = b_rvm(1,mu(j),kappa(j));
    end

% Parameter estimation for the bootstrap sample (Step 4.)
    if k == 1
        [ftm_star, mu_star, mvl_star] = mvlmn(theta_star,'rad');
        kappa_star = A1inv(mvl_star);
        p_star = [];
    else
        [param,err] = watsonkfit(theta_star,k);     % parameter estimation
        mu_star = param{1};
        kappa_star = param{2};
        p_star = param{3};
    end
    p_star(end+1) = 1 - sum(p_star);

    F_star = cell(1,k);    % mixed von Mises cdf
    cF_star = zeros(1,n);
    for t = 1:k
        F_star{t} = vmcdf(theta_star,mu_star(t),kappa_star(t),'rad');
        cF_star = cF_star + p_star(t) * F_star{t};
    end

    z = sort(cF_star,'ascend');  % goodness-of-fit
    T0 = (z - (2 * (1:n) - 1) / (2 * n)) .^ 2;
    T_star(bts) = sum(T0) - n * (mean(z) - 0.5) ^ 2 + 1 / (12 * n);
    save([resdir 'Tstar_' num2str(k) '.mat'],'T_star')
end

% Test k modes against more than k modes
pvalue = (1 / B) * sum(T<=T_star);