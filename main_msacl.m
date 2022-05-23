function main_msacl
%MS_MOD_TSC_MAIN Main wrapper for the examples presented in the Király &
% Hangya 'Model selection and clustering minefields review.

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   23-May-2022


%% Figure 1
% Examples of model selection problems in neuroscience
%--------------------------------------------------------------------------
% Panel a-b
% Using MAICE to choose between competing tuning curve gain models 
% Panel a left 
gain_factor = 2; % multiplicative factor applied for data simulation
additive_factor = 0; % additive factor applied for data simulation
simnum = 1; % number of simulations
tuning_model_sim(gain_factor,additive_factor,simnum);
% Panel a right
gain_factor = 2;
additive_factor = 0;
simnum = 100;
tuning_model_sim(gain_factor,additive_factor,simnum);
% Panel b
gain_factor = 1.5;
additive_factor = 1;
simnum = 1;
tuning_model_sim(gain_factor,additive_factor,simnum);

% Panel c
% Using BIC to choose the best fitting ARMA model
ARMA_simulation;

%Panel d-f
% Using information criteria and parametric bootstrap technique for choosing
% the number of modes in a distribution.
% Panel d
N = 250; % Number of neurons
D = 1; % Distance paramter (between the two closer modes) in radian
components =  [1,2,3,4]; % number of modes for the fitted distributions
B_simnum = 100; % number of bootstrap simulations 
warning('off','all')
phaselock_modes_simulation_param(N,D,components,B_simnum);
% Panel e
N = 250;
D = 0.85;
components =  [1,2,3,4];
B_simnum = 100;
phaselock_modes_simulation_param(N,D,components,B_simnum);
% Panel f left
N = 250;
D = 1:0.2:2;
components = [2,3];
dAIC = zeros(1,length(D)); % preallocate memory
dBIC = zeros(1,length(D));
for i = 1:length(D)
    [dAIC(i),dBIC(i)] = phaselock_modes_simulation_param(N,D(i),components,0);
end
figure
scatter (D,flip(dAIC),'g') 
hold on
scatter (D,flip(dBIC),'b')
line([2,3],[0,0])
% Panel f right
N = [250,2750,5250,7750];
D = 1.15;
components = [2,3];
dAIC = zeros(1,length(N)); % preallocate memory
dBIC = zeros(1,length(N));
for i = 1:length(N)
    [dAIC(i),dBIC(i)] = phaselock_modes_simulation_param(N(i),D,components,0);
end
figure
scatter (N,dAIC,'g')
hold on
scatter (N,dBIC,'b')
line([1,8000],[0,0])

%% Figure 2
% Examples of clustering problems in neuroscience
%--------------------------------------------------------------------------
% Fig 2c
% Hierarchical clustering of simulated peri-event time histogram data
psth_clust

% Fig 2d
% Clustering of human cells based on RNA-sequencing data
% load publicaly available data
% https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq
load('L:\Balint\NatNeuro_prespective\RNA_all.mat')
load('L:\Balint\NatNeuro_prespective\RNA.mat') % load marker genes
clusternum = 3; % required cluster nuber
PCA_num = 20; % Number of pricinpal components to be used
rnaseq_clustering(RNA_all,RNA,clusternum,PCA_num)

% Fig 2e 
% Spike sorting of simulated action potentials using k-means clustering.
rng(1) % initialize seed for random number generator
start_coo = [20,20,85,45;20,50,20,20]'; % starting centroid coodinates #1
kmeans_spikesorting_simulation(start_coo)
start_coo = [20,20,35,50;20,50,20,20]'; % starting centroid coodinates #2
kmeans_spikesorting_simulation(start_coo)
start_coo = [20,20,20,50;20,50,25,20]'; % starting centroid coodinates #3
kmeans_spikesorting_simulation(start_coo)

