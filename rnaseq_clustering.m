function rnaseq_clustering(expression,marker_expression,clusternum,PCA_num)
%RNASEQ_CLUSTERING    Clustering of cells based on gene expression
%RNASEQ_CLUSTERING(EXPRESSION,M_EXPRESSION,CLUSTERNUM,PCA_NUM) %  Performs
% agglomerative hierarhical clustering on RNA-sequencing data (EXPRESSION) 
% based on the first PCA_NUM principcal components to CLUSTERNUM clusters.
% Soft k-means clustering based on the gene expression for two clusters 
% (excitatory, inhibitory) is also performed for comparison.
% The expresson of marker genes (MARKER_EXPRESSION) is used for
% visualization. 
% EXPRESSION is a n-by-m matrix where n is the number of cells and m is the
% number of genes.
% MARKER_EXPRESSION is a n-by-m matrix where n is the number of cells and m
% is the number of marker genes.

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   29-Jun-2021


% The first PCA_num principal compenents of the gene expression is
% calculated
[~,PCA_all]=pca(expression(:,:));
PCA=PCA_all(:,1:PCA_num);

% The dendrogram of cells is prepared based on the principcal components
% and cutted at clusternum clusters
Dend = linkage(PCA,'average','euclidean'); 
cutoff = Dend(end-3+2,clusternum);

% The dendrogram is plotted with the clusters
subplot(3,1,1)
[~,~,leafOrder]=dendrogram(Dend,0,'Orientation','top','ColorThreshold',cutoff);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'xtick',[])
set(gca,'ytick',[])

% marker gene expression is plotted
subplot(3,1,2)
imagesc(marker_expression(leafOrder,:)')
set(gca,'xtick',[])
set(gca,'ytick',[])

% soft k-means clustering is performed and probabilities are plotted
[~,U] = fcm(expression(:,:),2);
subplot(3,1,3)
bar(U(:,leafOrder)', 'stacked')
ylim([0,1])
set(gca,'xtick',[])
set(gca,'ytick',[0,0.5,1])