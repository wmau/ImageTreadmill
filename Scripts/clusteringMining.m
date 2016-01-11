clear;
close all;

load('A0.01.mat');

animal = 'GCamp6f_45_treadmill';
date = '11_20_2015';
session = 2;

centroids = getNeuronCentroids(animal,date,session);

A(A>1) = 1; 
Au(Au>1) = 1; 
Au(eye(size(Au))~=0) = 0; 

v = GCModulMax3(A);

% [modules,inmodule] = louvain_community_finding(Au);
% v = cell2mat(inmodule);

tbl = tabulate(v); 
clusters = tbl(tbl(:,2)>2,1);
nClusters = length(clusters); 
inCluster = find(ismember(v,clusters));

null = pdist(centroids,'euclidean');

d = cell(nClusters,1);
H = zeros(nClusters,1); p = zeros(nClusters,1);
figure(1); scatter(centroids(:,1),centroids(:,2),'.'); hold on;
figure(2); hold on; 
for i=1:nClusters
    thisCluster = clusters(i);
    good = find(v==thisCluster);
    figure(1);
    scatter(centroids(good,1),centroids(good,2),1000,'.');
    
    d{i} = pdist(centroids(v==thisCluster,:),'euclidean'); 
    
    figure(2);
    ecdf(d{i});
    [H(i),p(i)] = kstest2(d{i},null);
    
    if H(i)
        h = get(gca,'children');
        set(h(1),'linestyle','--');
        disp(num2str(clusters(i)));
    end
end

ecdf(null);
h = get(gca,'children');
set(h(1),'linewidth',3);
