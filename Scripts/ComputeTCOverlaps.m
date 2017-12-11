clear;
loadMD;

fulldataset = [MD(292:299) MD(300:303) MD(305:308)];

animals = unique({fulldataset.Animal});
nAnimals = length(animals); 
B = 100;
reference = 'day1MappedTCs';

p = nan(nAnimals,3);
r = nan(nAnimals,3,B); 
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    j=1;
    for s=2:length(ssns)
        [~,~,p(a,j)] = TCEnsembleOverlap(fulldataset(ssns(1)),fulldataset(ssns(s)),...
            reference,false);
        
        for i=1:B
            [~,~,r(a,j,i)] = TCEnsembleOverlap(fulldataset(ssns(1)),fulldataset(ssns(s)),...
                reference,true);
        end
        j=j+1;
    end
end

figure('Position',[790 325 265 465]); hold on;
plot(p','color',[.7 .7 .7]); 
errorbar(mean(p),std(p)/2,'k','linewidth',3);
errorbar(mean(mean(r,3)),std(mean(r,3))/2,'linewidth',3,'color','r');
xlim([0.5 3.5]);
xlabel('Days later');
ylabel('Overlapping ensemble (% membership)');
set(gca,'tickdir','out');

nComparisons = size(p,2);
randomOverlap = zeros(B*nAnimals,nComparisons);
for i=1:nComparisons
    temp = squeeze(r(:,i,:));
    randomOverlap(:,i) = temp(:);
end

dayLag = [ones(nAnimals,1),2.*ones(nAnimals,1),3.*ones(nAnimals,1)];
dayLagShuffle = [ones(nAnimals*B,1),2.*ones(nAnimals*B,1),3.*ones(nAnimals*B,1)];
emp = ones(size(p));
shuffle = zeros(size(r));
[pval,tab,stats] = anovan([p(:); r(:)],{[emp(:);shuffle(:)],[dayLag(:);dayLagShuffle(:)]});