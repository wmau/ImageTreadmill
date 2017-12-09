clear
loadMD;

    %MD(292:295) = G45.
    %MD(296:299) = G48.
    %MD(300:304) = Bellatrix.
    %MD(305:309) = Polaris.
    
fulldataset = MD(292:309);
animals = unique({fulldataset.Animal});
nAnimals = length(animals);
colors = parula(nAnimals);

statType = 'si';
switch statType
    case 'ti', col = [0 .5 .5]; infoType = 'Temporal';
    case 'si', col = [.58 .44 .86]; infoType = 'Spatial';
end
PCcrit = .01;

saveBool = false;
folder = 'C:\Users\William Mau\Documents\Projects\Time Cell Imaging Summer 2015 -\Paper\Figures';
scattername = fullfile(folder,['Stable ',statType,' Scatter']);
linename = fullfile(folder,['Stable ',statType,' Drop']);

if saveBool
    c = input('Saving set to true. Are you sure you want to continue? (y/n)','s');
    
    if ~strcmp(c,'y')
        saveBool = false;
    end
end

nSessions = zeros(1,nAnimals);
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    nSessions(a) = length(ssns)-1;
end

[S,US] = deal(cell(1,3));
[SM,USM] = deal(cell(nAnimals,1));
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    nSessions = length(ssns);
    [SM{a},USM{a}] = deal(cell(nSessions,1));
    
    for s=2:4
        cd(fulldataset(ssns(1)).Location);
        load('FinalOutput.mat','NumNeurons');
        
        switch statType
            case 'ti',neurons = getTimeCells(fulldataset(ssns(1))); 
            case 'si',neurons = getPlaceCells(fulldataset(ssns(1)),PCcrit);
        end
        
        neurons = EliminateUncertainMatches([fulldataset(ssns(1)) fulldataset(ssns(s))],...
            neurons);
        
        switch statType
            case 'si', corrs = CorrPlaceFields(fulldataset(ssns(1)),fulldataset(ssns(s)),neurons);
            case 'ti', corrs = CorrTrdmllTrace(fulldataset(ssns(1)),fulldataset(ssns(s)),neurons);
        end
        
        stblcrit = .01/length(neurons);
        stable = intersect(find(corrs(:,2) < stblcrit),neurons);
        unstable = intersect(find(corrs(:,2) > stblcrit | isnan(corrs(:,2))),neurons);
        
        temp = msStats(fulldataset([ssns(1) ssns(s)]),statType,1:NumNeurons);
        S{s-1} = [S{s-1}; temp(stable,:)];
        US{s-1} = [US{s-1}; temp(unstable,:)];
        
        if length(stable) > 3
            sTemp = temp(stable,:); 
            sTemp(sTemp(:,2)==0,:) = [];
            SM{a}{s} = mean(sTemp,1);
        else, SM{a}{s} = nan(1,2);
        end
        
        if length(unstable) > 3
            uTemp = temp(unstable,:);
            uTemp(uTemp(:,2)==0,:) = [];
            USM{a}{s} = mean(uTemp,1);
        else, USM{a}{s} = nan(1,2); 
        end
        
        
    end
    
end

d = length(S); 
c = jet(d); 
figure; hold on;
for i=d:-1:1
    s(i) = scatter(S{i}(:,1),S{i}(:,2),[],c(i,:),'filled');
end
alpha(s,.5); 
l = lsline; 

for i=1:d
    l(i).Color = c(i,:);
end

% for i=1:d
%     diffs{i} = diff(S{i},[],2);
%     grps{i} = i*ones(1,length(diffs{i}));
% end
% D = cell2mat(diffs');
% grps = cell2mat(grps);
% 
% figure;
% [~,~,stats] = anovan(D',grps');
% multcompare(stats);
% 
% for i=1:d
%     M(i) = mean(D(grps==i));
%     SEM(i) = std(D(grps==i))/sqrt(length(D(grps==i)));
% end
% figure;
% errorbar(M,SEM,'linewidth',3,'color',col);
% 
% xlim([.5 3.5]);
% set(gca,'xtick',[1:1:3],'fontsize',15,'linewidth',4);
% xlabel('Days from Reference'); 
% ylabel(['Delta ',infoType,' Info. [bits]']);

%% Look only at quartiles. 
for i=1:d
    allInfos = [S{i}; US{i}];
    upperQuartile = quantile(allInfos(:,1),.75);
    lowerQuartile = quantile(allInfos(:,1),.25);
    in = allInfos(:,1)>upperQuartile;
    %in = allInfos(:,1)<lowerQuartile;
    diffs{i} = diff(allInfos(in,:),[],2);
    grps{i} = i*ones(1,length(diffs{i}));
end
D = cell2mat(diffs');
grps = cell2mat(grps);

figure;
[~,~,stats] = anovan(D',grps');
multcompare(stats);

for i=1:d
    M(i) = mean(D(grps==i));
    SEM(i) = std(D(grps==i))/sqrt(length(D(grps==i)));
end
figure;
errorbar(M,SEM,'linewidth',3,'color',col);

xlim([.5 3.5]);
set(gca,'xtick',[1:1:3],'fontsize',15,'linewidth',4);
xlabel('Days from Reference'); 
ylabel(['Delta ',infoType,' Info. [bits]']);
