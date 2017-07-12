%% Description.
% Determine whether entering/exiting cells have more or less displacement
% than their stable peers. 
%

%% Set up.
clear; 
loadMD;

folder = 'C:\Users\William Mau\Documents\Projects\Time Cell Imaging Summer 2015 -\Paper\Figures\Supplementals';
fn = fullfile(folder,'Stability Drifts');
cellType = 'time';
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

animals = unique({mds.Animal});
nAnimals = length(animals);

%% Loop through animals and sessions. 
[stableCentroidDrifts,stableOrientationDrifts,...
    exitCentroidDrifts,exitOrientationDrifts,...
    enterCentroidDrifts,enterOrientationDrifts] = deal(cell(nAnimals,1));
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

    sRows = []; exRows = []; entRows = [];
    [stableCentroidDrifts{a},stableOrientationDrifts{a},...
        exitCentroidDrifts{a},exitOrientationDrifts{a},...
        enterCentroidDrifts{a},enterOrientationDrifts{a}] = deal(cell(1,length(ssns)-1));
    for s = 1:length(ssns)-1
        cd(mds(ssns(s)).Location);
        %stable = union(Sstable{a}{s},Tstable{a}{s});
        %unstable = union(Sunstable{a}{s},Tunstable{a}{s});

        %Get cell stability status. 
        [stable,exiting,entering] = CellStabilityStatus(mds(ssns(s)),mds(ssns(s+1)),cellType);

        %Get rows in the map matrix that correspond to each cell category.
        [~,mapRows,mapCols(s)] = msMatchCells(mapMD,mds(ssns(s)),stable,false);
        sRows = [sRows; mapRows];
        [~,mapRows,mapCols(s)] = msMatchCells(mapMD,mds(ssns(s)),exiting,false);
        exRows = [exRows; mapRows];
        [~,mapRows,mapCols(s)] = msMatchCells(mapMD,mds(ssns(s)),entering,false);
        entRows = [entRows; mapRows];
            
    end
    
    %Get cell numbers for each category. 
    sMAP{a} = MAP{a}(sRows,mapCols);
    exMAP{a} = MAP{a}(exRows,mapCols);
    entMAP{a} = MAP{a}(entRows,mapCols);
    
    %Get cell distance for each cell category. 
    for s = 1:length(ssns)-1
        %stable = union(Sstable{a}{s},Tstable{a}{s});
        %unstable = union(Sunstable{a}{s},Tunstable{a}{s});
        
        %Get cell category. 
        [stable,exiting,entering] = CellStabilityStatus(mds(ssns(s)),mds(ssns(s+1)),cellType);
        
        %Get cell distances. 
        %Stable cells.
        reg_stats = neuron_reg_qc(mds(ssns(s)),...
            mds(ssns(s+1)),'neurons',intersect(sMAP{a}(:,s),stable));
        stableCentroidDrifts{a}{s} = reg_stats.cent_d;
        stableM(i) = mean(reg_stats.cent_d);
        stableOrientationDrifts{a}{s} =  reg_stats.orient_diff;

        %Exiting cells.
        reg_stats = neuron_reg_qc(mds(ssns(s)),...
            mds(ssns(s+1)),'neurons',intersect(exMAP{a}(:,s),exiting));
        exitCentroidDrifts{a}{s} = reg_stats.cent_d;
        exitM(i) = mean(reg_stats.cent_d);
        exitOrientationDrifts{a}{s} = reg_stats.orient_diff;

        %Entering cells. 
        reg_stats = neuron_reg_qc(mds(ssns(s)),...
            mds(ssns(s+1)),'neurons',intersect(entMAP{a}(:,s),entering));
        enterCentroidDrifts{a}{s} = reg_stats.cent_d;
        enterM(i) = mean(reg_stats.cent_d);
        enterOrientationDrifts{a}{s} = reg_stats.orient_diff;

        i=i+1;
                   
    end
end

%% Plot. 
i=1;
figure('Position',[130 240 250 440]); hold on;
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{mds.Animal}));
    
    %Individual sessions.
    for s=1:length(ssns)-1
        plot([1,2,3],[stableM(i),exitM(i),enterM(i)],'-','color',[.7 .7 .7]);
        
        i=i+1;
    end
end
%Mean and SEM. 
m = [nanmean(stableM), nanmean(exitM), nanmean(enterM)];
sem = [nanstd(stableM)/sqrt(length(stableM)) nanstd(exitM)/sqrt(length(exitM)) ...
    nanstd(enterM)/sqrt(length(enterM))];
errorbar([1,2,3],m,sem,'color','k','linewidth',5);
set(gca,'tickdir','out','linewidth',4,'fontsize',12,'xlim',[0.5 3.5],...
    'xtick',[1:3],'xticklabel',{'Stable','Exiting','Entering'});
ylabel('Mean ROI displacement (m)','fontsize',15);

if saveBool
    print(fn,'-dpdf');
end

%% Histograms. 
s = []; ex = []; ent = [];
for i=1:4
    s = [s; cell2mat(stableCentroidDrifts{i}')];
    ex = [ex; cell2mat(exitCentroidDrifts{i}')];
    ent = [ent; cell2mat(enterCentroidDrifts{i}')];
end
figure; hold on

[~,p] = kstest2(s,ex);
[N,edges] = histcounts(s,'binwidth',.1);
[N,edges] = hist(s,edges,'normalization','probability');
N = N./sum(N);
stairs(edges,N,'linewidth',4);

[N,edges] = histcounts(ex,'binwidth',.1);
[N,edges] = hist(ex,edges,'normalization','probability');
N = N./sum(N);
stairs(edges,N,'linewidth',4);

[N,edges] = histcounts(ent,'binwidth',.1);
[N,edges] = hist(ent,edges,'normalization','probability');
N = N./sum(N);
stairs(edges,N,'linewidth',4);

set(gca,'tickdir','out','linewidth',4,'fontsize',15);
xlim([0 3.3])
legend({'Stable','Exiting','Entering'});
ylabel('Proportion');
xlabel('Centroid drift [microns]');
title(['KS P = ',num2str(p)]); 

d = [stableM exitM enterM];
grps = [zeros(1,14) ones(1,14) 2*ones(1,14)];
[~,~,stats] = anovan(d,{grps})