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
mds = MD(292:309);      

[~,~,Sstable,Sunstable] = PartitionStats(mds,'place','si');
[~,~,Tstable,Tunstable] = PartitionStats(mds,'time','ti');

animals = unique({mds.Animal});
nAnimals = length(animals);
c = parula(nAnimals);

[stableCentroidDrifts,stableOrientationDrifts,...
    unstableCentroidDrifts,unstableOrientationDrifts,...
    stableOverlap,unstableOverlap] = deal(cell(nAnimals,1));
i=1;
for a=1:nAnimals
    %Get all the sessions for this animal.
    ssns = find(strcmp(animals{a},{mds.Animal}));

    %Load registration matrix.
    mapMD = getMapMD(mds(ssns));
    cd(mapMD.Location);
    load('batch_session_map.mat');
    MAP{a} = batch_session_map.map(:,2:end);

    mapCols = zeros(length(ssns)-1,1);
        
    sRows = []; usRows = [];
    [stableCentroidDrifts{a},stableOrientationDrifts{a},...
        unstableCentroidDrifts{a},unstableOrientationDrifts{a},...
        stableOverlap{a},unstableOverlap{a}] = deal(cell(1,length(ssns)-1));
    for s = 1:length(ssns)-1
        cd(mds(ssns(s)).Location);
        %stable = union(Sstable{a}{s},Tstable{a}{s});
        %unstable = union(Sunstable{a}{s},Tunstable{a}{s});

        stable = Sstable{a}{s};
        unstable = Sunstable{a}{s};

        [~,mapRows,mapCols(s)] = msMatchCells(mds(ssns(s)),stable,false);
        sRows = [sRows; mapRows];
        [~,mapRows,mapCols(s)] = msMatchCells(mds(ssns(s)),unstable,false);
        usRows = [usRows; mapRows];
            
    end
    
    sMAP{a} = MAP{a}(sRows,mapCols);
    usMAP{a} = MAP{a}(usRows,mapCols);
    
    for s = 1:length(ssns)-1
        %stable = union(Sstable{a}{s},Tstable{a}{s});
        %unstable = union(Sunstable{a}{s},Tunstable{a}{s});
        
        stable = Sstable{a}{s};
        unstable = Sunstable{a}{s};

        reg_stats = neuron_reg_qc(mds(ssns(s)),...
            mds(ssns(s+1)),'neurons',intersect(sMAP{a}(:,s),stable));
        stableCentroidDrifts{a}{s} = reg_stats.cent_d;
        stableM(i) = mean(reg_stats.cent_d);
        stableOrientationDrifts{a}{s} =  reg_stats.orient_diff;
        stableOverlap{a}{s} = reg_stats.overlap;

        reg_stats = neuron_reg_qc(mds(ssns(s)),...
            mds(ssns(s+1)),'neurons',intersect(usMAP{a}(:,s),unstable));
        unstableCentroidDrifts{a}{s} = reg_stats.cent_d;
        unstableM(i) = mean(reg_stats.cent_d);
        unstableOrientationDrifts{a}{s} = reg_stats.orient_diff;
        unstableOverlap{a}{s} = reg_stats.overlap;


        i=i+1;
                   
    end
end

i=1;
figure('Position',[130 240 250 440]); hold on;
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{mds.Animal}));
    
    for s=1:length(ssns)-1
        plot([1,2],[stableM(i),unstableM(i)],'-','color',c(a,:));
        
        i=i+1;
    end
end
m = [mean(stableM), mean(unstableM)];
sem = [std(stableM)/sqrt(length(stableM)) std(unstableM)/sqrt(length(unstableM))];
errorbar([1,2],m,sem,'color','k','linewidth',5);

[pCD,~,z] = signrank(stableM,unstableM,'method','approximate');
title(['Z = ',num2str(z.zval),' P = ',num2str(pCD)]);
set(gca,'xticklabel',{'Stable','Unstable'},'tickdir','out',...
    'linewidth',4,'xtick',[1:2],'fontsize',15);
xlim([0.5,2.5]); ylim([0,3]);
ylabel('Mean Centroid Drift [microns]');

if saveBool
    print(fn,'-dpdf');
end

s = []; us = [];
for i=1:4
    s = [s; cell2mat(stableCentroidDrifts{i}')];
    us = [us; cell2mat(unstableCentroidDrifts{i}')];
end
[~,p] = kstest2(s,us);
[N,edges] = histcounts(s,'binwidth',.1);
[N,edges] = hist(s,edges,'normalization','probability');
N = N./sum(N);
figure;
stairs(edges,N,'linewidth',4);
hold on
[N,edges] = histcounts(us,'binwidth',.1);
[N,edges] = hist(us,edges,'normalization','probability');
N = N./sum(N);
stairs(edges,N,'linewidth',4);
set(gca,'tickdir','out','linewidth',4,'fontsize',15);
xlim([0 3.3])
legend({'Stable','Unstable'});
ylabel('Proportion');
xlabel('Centroid drift [microns]');
title(['KS P = ',num2str(p)]); 