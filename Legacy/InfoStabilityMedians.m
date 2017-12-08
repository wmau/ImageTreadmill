%% Set up.
    clear;
    
    %Load MD entries. 
    loadMD;
    %MD(292:295) = G45.
    %MD(296:299) = G48.
    %MD(300:304) = Bellatrix.
    %MD(305:309) = Polaris.
    fulldataset = MD(292:309);      

    %Some initial variables. 
    teal = [0 .5 .5];
    purple = [.58 .44 .86];
    
%% Save information
    saveBool = true;
    folder = 'C:\Users\William Mau\Documents\Projects\Time Cell Imaging Summer 2015 -\Paper\Figures\Supplementals\Alternative Info Stability Analyses';
    TimeTI = fullfile(folder,'Stable Time TI Medians');
    TimeSI = fullfile(folder,'Stable Time SI Medians');
    TimeFR = fullfile(folder,'Stable Time FR Medians');
    PlaceSI = fullfile(folder,'Stable Place SI Medians');
    PlaceTI = fullfile(folder,'Stable Place TI Medians');
    PlaceFR = fullfile(folder,'Stable Place FR Medians'); 
    
    if saveBool
        c = input('Saving set to true. Are you sure you want to continue? (y/n)','s');
    
        if ~strcmp(c,'y')
            saveBool = false;
        end
    end
    
    ParseInfoStability(fulldataset,'time','ti');
    ylabel('Median Temp. Info. [z-scored bits]');
    if saveBool
        print(TimeTI,'-dpdf');
    end
    
    ParseInfoStability(fulldataset,'time','si');
    ylabel('Median Spat. Info. [z-scored bits]');
    if saveBool
        print(TimeSI,'-dpdf');
    end
    
    ParseInfoStability(fulldataset,'time','fr');
    ylabel('Normalized Ca Event Frequency');
    if saveBool
        print(TimeFR,'-dpdf');
    end
    
    ParseInfoStability(fulldataset,'place','si');
    ylabel('Median Spat. Info. [z-scored bits]');
    if saveBool
        print(PlaceSI,'-dpdf');
    end
    
    ParseInfoStability(fulldataset,'place','ti');
    ylabel('Median Temp. Info. [z-scored bits]');
    if saveBool
        print(PlaceTI,'-dpdf');
    end
    
    ParseInfoStability(fulldataset,'place','fr');
    ylabel('Normalized Ca Event Frequency');
    if saveBool
        print(PlaceFR,'-dpdf');
    end
    
function ParseInfoStability(fulldataset,stabilityType,infoType)
    animals = unique({fulldataset.Animal});
    nAnimals = length(animals);
    c = parula(nAnimals);

    [stats,~,stable,unstable] = PartitionStats(fulldataset,stabilityType,infoType);
    for a=1:nAnimals
        ssns = find(strcmp(animals{a},{fulldataset.Animal}));
        
        [i,j] = deal(1);
        for s=1:length(ssns)-1
            sExtent = length(stable{a}{s}); 
            usExtent = length(unstable{a}{s});
            sM{a}(s) = median(stats.stable{a}(i:i+sExtent-1));
            usM{a}(s) = median(stats.unstable{a}(j:j+usExtent-1)); 
            
            i = i+sExtent;
            j = j+usExtent;
        end
    end
    
    sM = cell2mat(sM);
    usM = cell2mat(usM);
    p = ranksum(sM,usM);
    i = 1;
    figure('Position',[790 220 250 440]);
    hold on
    for a=1:nAnimals
        ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    
        for s=1:length(ssns)-1
            hold on;
            plot([1,2],[sM(i),usM(i)],'-o','color',c(a,:),'linewidth',3);

            i=i+1;
        end
    end
    set(gca,'xticklabel',{'Stable','Unstable'});
    set(gca,'xtick',[1:2],'tickdir','out');
    title(['P = ',num2str(p)]);
    xlim([0.5,2.5]);
end
    