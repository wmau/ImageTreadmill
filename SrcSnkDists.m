function [d,dnull] = SrcSnkDists(md,A)
%[d,dnull] = SrcSnkDists(md,A)
%
%   

%% 
    cd(md.Location); 
    B = 1000;
    
    %Neuron centroids. 
    cntrds = getNeuronCentroids(md);
    
    %List of all neurons active on the treadmill.
    [nNrns,trdmllNrns] = nNeuronsActiveonTM(md); 
    
    %Get edge list.
    el = adj2edgeL(A); el = el(:,1:2);
    nE = size(el,1); 
    
    %Number of sink / target neurons. 
    snks = unique(el(:,2));
    nSnks = length(snks);
    
    %Get distances to random cells.
    dnull = zeros(nSnks,B);
    for s = 1:nSnks
        smpl = randsample(trdmllNrns(trdmllNrns~=snks(s)),B,true);
        dnull(s,:) = fastdist(cntrds(smpl,:),cntrds(snks(s),:));
    end
    
    %Get real distances. 
    d = nan(size(A));
    p = nan(size(A));
    for e = 1:nE 
        src = el(e,1);
        snk = el(e,2); 
        
        %Pairwise distance square matrix. And associated p-values. 
        d(src,snk) = fastdist(cntrds(src,:),cntrds(snk,:));
        p(src,snk) = sum(d(src,snk) > dnull(snks==snk,:))/B;
    end
    
    %Plot distance histogram.
    figure; hold on ;
    histogram(d(~isnan(d)),'binwidth',20,'normalization','probability','facecolor','y');
    histogram(dnull(:),'binwidth',20,'normalization','probability','facecolor','c');
    xlabel('Distances [pixels]'); 
    ylabel('Proportion'); 
    [~,pval] = kstest2(d(~isnan(d)),dnull(:)); 
    title(['KS p = ',num2str(pval)]);
    
end

function d = fastdist(comp,ref)
%
%
%
    
%%  
    nComps = size(comp,1);
    d = nan(1,nComps);
    
    for c = 1:nComps
       d(c) =  sqrt(((ref(1) - comp(c,1)) * (ref(1) - comp(c,1))) + ...
           ((ref(2) - comp(c,2)) * (ref(2) - comp(c,2))));
    end
end
