%% Description
%
%   Examines peaks of new time cells the day before and correlates them. If
%   the correlation is significant, it suggests that the peaks were in the
%   same place. 
%

%% Set up.
    clear;
    
    %Load MD entries. 
    loadMD;
    %MD(292:295) = G45.
    %MD(296:299) = G48.
    %MD(300:304) = Bellatrix.
    %MD(305:309) = Polaris.
    fulldataset = MD(292:309);      

    S1 = [];
    S2 = [];
    
    purple = [.58 .44 .86];
    
    animals = unique({fulldataset.Animal});
    nAnimals = length(animals);
    
    for a=1:nAnimals
        ssns = find(strcmp(animals{a},{fulldataset.Animal})); 
        nSessions = length(ssns)-1;
        
        for s=1:nSessions
            [a,b] = RankNewPCs(fulldataset(ssns(s)),fulldataset(ssns(s+1)),'plotit',false);
            
            S1 = [S1; a];
            S2 = [S2; b]; 
        end
    end
    
    [R,p] = corr(S1,S2,'type','spearman');
    figure;
    scatter(S1,S2,20,purple,'filled');
    title(['R = ',num2str(R),' p = ',num2str(p)]);
    line([0:max([S1; S2])],[0:max([S1; S2])],'linestyle','-.','color','k');
    xlabel('Day 1 Peak'); 
    ylabel('Day 2 Peak');
    set(gca,'linewidth',2,'tickdir','out');