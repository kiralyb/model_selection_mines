function kmeans_spikesorting_simulation(start)
%KMEANS_SPIKESORTING_SIMULATION     Performs k-means spike sorting on
%simulated data from starting centroid locations defined in START and
%applys the Elbow method to find the optimal number of clusters

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   29-Jun-2021

% simulate x-y features of action potentials 
% predefine cluster centers
x = [ones(1,100) * 20,ones(1,25) * 38,ones(1,25) * 53,ones(1,30) * 20,ones(1,40) * 20,80];
y = [ones(1,100) * 20,ones(1,50) * 20,ones(1,30) * 30,ones(1,40) * 50,20];
Clustnum = 4; % number of clusters
max_Clustnum = 10; % maximal number of clusters tested

% add gaussian variance from the centers
%rng(1)
for i = 1:length(x) 
    x(i) = x(i) + normrnd(0,5);
    y(i) = y(i) + normrnd(0,5);
end

% perform k-means clustering with 100 repetitions for each k in the range 
% of 1 to 10, and measure the inter-cluster squared Euclidean 
% distance for each cluster 
D = zeros(max_Clustnum);
for j = 1:max_Clustnum
    [~,~,dist] = kmeans([x;y]',j,'Replicates',100);
    D(1:length(dist),j) = dist';
end

% plot the summed inter-cluster squared Euclidean distance for each k to
% allow visual identification of the 'elbow'
figure
plot(sum(D));

% plot the result of k-means with the given initial centroid locations
[c,~] = kmeans([x;y]',Clustnum,'start',start);
figure
scatter(x,y,[],c,'x')
hold on
plot(start(:,1),start(:,2),'kx','MarkerSize',15,'LineWidth',3) 
