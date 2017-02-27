clear
close all;
loadMD; 

fulldataset = MD(292:309);

animals = unique({fulldataset.Animal});
nAnimals = length(animals); 

d = [];
o = [];
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    
    for s=1:length(ssns)-1
        s1 = fulldataset(ssns(s));
        s2 = fulldataset(ssns(s+1));
        
        regstats = neuron_reg_qc(s1,s2);
        d = [d; regstats.cent_d];
        o = [o; regstats.orient_diff];
    end
end

figure;
histogram(d,40,'normalization','probability','edgecolor','none');
xlabel('Centroid distance [microns]');
ylabel('Proportion');
set(gca,'tickdir','out');

figure;
histogram(o,40,'normalization','probability','edgecolor','none');
xlabel('Orientation difference [degrees]');
ylabel('Proportion');
set(gca,'tickdir','out');

