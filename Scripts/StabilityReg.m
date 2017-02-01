clear; 
loadMD;

folder = 'C:\Users\William Mau\Documents\Projects\Time Cell Imaging Summer 2015 -\Paper\Figures\Supplementals';
fn = fullfile(folder,'Stability Drifts');

saveBool = false;
if saveBool
    c = input('Saving set to true. Are you sure you want to continue? (y/n)','s');

    if ~strcmp(c,'y')
        saveBool = false;
    end
end

%MD(292:295) = G45.
%MD(296:299) = G48.
%MD(300:304) = Bellatrix.
%MD(305:309) = Polaris.
fulldataset = MD(292:309);      

[~,~,Sstable,Sunstable] = PartitionStats(fulldataset,'place','si');
[~,~,Tstable,Tunstable] = PartitionStats(fulldataset,'time','ti');

animals = unique({fulldataset.Animal});
nAnimals = length(animals);
c = parula(nAnimals);

[stableCentroidDrifts,stableOrientationDrifts,...
    unstableCentroidDrifts,unstableOrientationDrifts,...
    stableOverlap,unstableOverlap] = deal(cell(nAnimals,1));
i=1;
figure;
for a=1:nAnimals
    %Get all the sessions for this animal.
    ssns = find(strcmp(animals{a},{fulldataset.Animal})); 
    
    [stableCentroidDrifts{a},stableOrientationDrifts{a},...
        unstableCentroidDrifts{a},unstableOrientationDrifts{a},...
        stableOverlap{a},unstableOverlap{a}] = deal(cell(1,length(ssns)-1));
    for s = 1:length(ssns)-1
            cd(fulldataset(ssns(s)).Location);
            stable = union(Sstable{a}{s},Tstable{a}{s});
            unstable = union(Sunstable{a}{s},Tunstable{a}{s});
                      
%             stable = Sstable{a}{s};
%             unstable = Sunstable{a}{s};
            
            reg_stats = neuron_reg_qc(fulldataset(ssns(s)),...
                fulldataset(ssns(s+1)),'neurons',stable);
            stableCentroidDrifts{a}{s} = reg_stats.cent_d;
            stableM(i) = median(reg_stats.cent_d);
            stableOrientationDrifts{a}{s} =  reg_stats.orient_diff;
            stableOverlap{a}{s} = reg_stats.overlap;
            
            reg_stats = neuron_reg_qc(fulldataset(ssns(s)),...
                fulldataset(ssns(s+1)),'neurons',unstable);
            unstableCentroidDrifts{a}{s} = reg_stats.cent_d;
            unstableM(i) = median(reg_stats.cent_d);
            unstableOrientationDrifts{a}{s} = reg_stats.orient_diff;
            unstableOverlap{a}{s} = reg_stats.overlap;
            
            
            i=i+1;
                   
    end
end

i=1;
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    
    for s=1:length(ssns)-1
        hold on;
        plot([1,2],[stableM(i),unstableM(i)],'-o','color',c(a,:),'linewidth',3);
        
        i=i+1;
    end
end

[pCD] = ranksum(stableM,unstableM);
title(['P = ',num2str(pCD)]);
% [~,pCD] = kstest2(stableCentroidDrifts,unstableCentroidDrifts);
% [~,pOD] = kstest2(stableOrientationDrifts,unstableOrientationDrifts);
% [~,pCorr] = ttest2(stableOverlap,unstableOverlap);
% 
% figure; 
% subplot(1,3,1); 
% ecdf(stableCentroidDrifts); 
% hold on;
% ecdf(unstableCentroidDrifts);
% title(['p = ',num2str(pCD)]);
% 
% subplot(1,3,2);
% ecdf(stableOrientationDrifts); 
% hold on;
% ecdf(unstableOrientationDrifts);
% title(['p = ',num2str(pOD)]);
% 
% subplot(1,3,3);
% ecdf(stableOverlap);
% hold on;
% ecdf(unstableOverlap);
% title(['p = ',num2str(pCorr)]);
set(gca,'xticklabel',{'Stable','Unstable'});
set(gca,'xtick',[1:2]);
xlim([0.5,2.5]); ylim([0,3]);
ylabel('Median Centroid Drift [microns]');

if saveBool
    print(fn,'-dpdf');
end