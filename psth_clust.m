function psth_clust
%PSTH_CLUST     Performing dimensionality reduction with principal component
% analysis and hierarchical clustering on simulated PSTH data. Clustered 
% PSTHs, the scatter plot of the PSTHs in the PCA space and the dendrogram is plotted. 

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   29-Jun-2021

% simualte data
y = data_setup(20,50,15);

% define parameters
PCA_num = 2; %number of principal comonents to use
clustnum = 4; %number of clusters

% perform clustering
[~,PCA1] = pca(y'); % running pricipal component analysis
PCA2 = PCA1(:,1:PCA_num); % extracting the first two pricinpal components 
Dend = linkage(PCA2,'complete','euclidean'); % calculating the dendrogram of psths
Clusters = cluster(Dend,'maxclust',clustnum); % cutting the tree at clusternum clusters
cutoff = Dend(end-clustnum+2,3); % cutting the tree at clusternum clusters
D = pdist(PCA2); % euclidean distrance between points in the artifical space
leafOrder2 = optimalleaforder(Dend,D); % optimal order for plotting

% plot results 
figure
subplot(1,2,1)
imagesc(y(:,leafOrder2)') %plotting sorted psth matrix 
subplot(1,2,2)
dendrogram(Dend,0,'Orientation','right','reorder',leafOrder2,'ColorThreshold',cutoff); %plotting dendrogram
set(gca,'Ydir','reverse');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
figure 
scatter(PCA2(:,1),PCA2(:,2),50,Clusters*100) %plotting the PCA space


function y = data_setup (D,M,S)
y = zeros(35,M);
y(15,1:D)=y(15,1:D)+3;
y(16,1:D)=y(16,1:D)+10;
y(17,1:D)=y(17,1:D)+7;
y(18,1:D)=y(18,1:D)+2;
y(25,1:D)=y(25,1:D)+7;
y(26,1:D)=y(26,1:D)+4;
y(27,1:D)=y(27,1:D)+1;

y(15,D+1)=y(15,D+1)+3;
y(16,D+1)=y(16,D+1)+10;
y(17,D+1)=y(19,D+1)+7;
y(18,D+1)=y(18,D+1)+5;
y(19,D+1)=y(19,D+1)+5;
y(20,D+1)=y(20,D+1)+5;
y(21,D+1)=y(21,D+1)+5;
y(22,D+1)=y(22,D+1)+5;
y(23,D+1)=y(23,D+1)+3;
y(24,D+1)=y(24,D+1)+3;
y(25,D+1)=y(25,D+1)+9;
y(26,D+1)=y(26,D+1)+4;
y(27,D+1)=y(27,D+1)+1;

y(15,D+2:D+S+1)=y(15,D+2:D+S+1)+1;
y(16,D+2:D+S+1)=y(16,D+2:D+S+1)+8;
y(17,D+2:D+S+1)=y(17,D+2:D+S+1)+8;
y(18,D+2:D+S+1)=y(18,D+2:D+S+1)+7;
y(19,D+2:D+S+1)=y(19,D+2:D+S+1)+6;
y(20,D+2:D+S+1)=y(20,D+2:D+S+1)+6;
y(21,D+2:D+S+1)=y(21,D+2:D+S+1)+6;
y(22,D+2:D+S+1)=y(22,D+2:D+S+1)+6;
y(23,D+2:D+S+1)=y(23,D+2:D+S+1)+5;
y(24,D+2:D+S+1)=y(24,D+2:D+S+1)+5;
y(25,D+2:D+S+1)=y(25,D+2:D+S+1)+5;
y(26,D+2:D+S+1)=y(26,D+2:D+S+1)+3;
y(27,D+2:D+S+1)=y(29,D+2:D+S+1)+3;
y(28,D+2:D+S+1)=y(28,D+2:D+S+1)+2;
y(29,D+2:D+S+1)=y(29,D+2:D+S+1)+1;
y(30,D+2:D+S+1)=y(30,D+2:D+S+1)+1;
y(31,D+2:D+S+1)=y(31,D+2:D+S+1)+1;

y(15,D+S+2:M-1)=y(15,D+S+2:M-1)-1;
y(16,D+S+2:M-1)=y(16,D+S+2:M-1)-4;
y(17,D+S+2:M-1)=y(17,D+S+2:M-1)-8;
y(18,D+S+2:M-1)=y(18,D+S+2:M-1)-7;
y(19,D+S+2:M-1)=y(19,D+S+2:M-1)-3;
y(20,D+S+2:M-1)=y(20,D+S+2:M-1)-2;
y(21,D+S+2:M-1)=y(21,D+S+2:M-1)-2;
y(22,D+S+2:M-1)=y(22,D+S+2:M-1)-1;

p=2.0;
y(15,M)=y(15,M)-1*p;
y(16,M)=y(16,M)-10*p;
y(17,M)=y(19,M)-9*p;
y(18,M)=y(18,M)-7*p;
y(19,M)=y(19,M)-3*p;
y(20,M)=y(20,M)-2*p;
y(21,M)=y(21,M)-1*p;

y=y./2;
rng(3);
y=awgn(y,1);
