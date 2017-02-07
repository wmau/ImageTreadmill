clear
loadMD
%MD(292:295) = G45.
%MD(296:299) = G48.
%MD(300:304) = Bellatrix.
%MD(305:309) = Polaris.
fulldataset = MD(292:309);   
animals = unique({fulldataset.Animal});
nAnimals = length(animals);

[new,r,N,R] = deal([]);
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    nSessions = length(ssns)-1;
    
    for s=1:nSessions
        [stat,randomSample] = NewTimeCellStats(fulldataset(ssns(s)),...
            fulldataset(ssns(s+1)),'ti');
        
        new = [new; stat];
        r = [r; randomSample];
        
        N = [N; median(new)];
        R = [R; median(r)];
    end
end

figure;
histogram(new,'normalization','probability','binwidth',.1,'edgecolor','none');
hold on
histogram(r,'normalization','probability','binwidth',.1,'edgecolor','none');

figure; hold on;
for i=1:length(N)
    plot([1 2],[N(i) R(i)],'o-','linewidth',3);
end
p = ranksum(N,R); 
title(['P = ',num2str(p)]);