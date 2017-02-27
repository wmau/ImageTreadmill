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
    
%% 
    [D,stabilityStatus] = StabilityLocations(fulldataset,'time');
    [s,r,S,R,c] = UnpackDs(D,stabilityStatus); 
    
%     figure; hold on;
%         histogram(s,'binwidth',.05,'edgecolor','none','normalization','probability');
%         histogram(r,'binwidth',.05,'edgecolor','none','normalization','probability');
%         [~,p] = kstest2(s,r); 
%             title(['P = ',num2str(p)]);
%             set(gca,'tickdir','out');
%             xlabel('Norm. Distance (z-microns)');
%             ylabel('Proportion');
            
    figure('Position',[790 220 250 440]); hold on;
        m = [mean(S) mean(R)];
        sem = [std(S)/sqrt(length(S)) std(R)/sqrt(length(R))];
        for i=1:length(S)
            plot([1,2],[S(i),R(i)],'-','color',c(i,:));
        end
        errorbar([1,2],m,sem,'linewidth',5,'color','k');
        p = signrank(S,R);
        set(gca,'xticklabel',{'Stable','Random'});
        set(gca,'xtick',[1:2],'tickdir','out');
        xlim([0.5,2.5]);
        ylabel('Mean Interneuronal Distance [microns]');
        title(['P = ',num2str(p)]);
        
%% 
    [D,stabilityStatus] = StabilityLocations(fulldataset,'place');
    [s,r,S,R,c] = UnpackDs(D,stabilityStatus); 
    
%     figure; hold on;
%         histogram(s,'binwidth',.05,'edgecolor','none','normalization','probability');
%         histogram(r,'binwidth',.05,'edgecolor','none','normalization','probability');
%         [~,p] = kstest2(s,r); 
%             title(['P = ',num2str(p)]);
%             set(gca,'tickdir','out');
%             xlabel('Norm. Distance (z-microns)');
%             ylabel('Proportion');
            
    figure('Position',[790 220 250 440]); hold on;
        m = [mean(S) mean(R)];
        sem = [std(S)/sqrt(length(S)) std(R)/sqrt(length(R))];
        for i=1:length(S)
            plot([1,2],[S(i),R(i)],'-','color',c(i,:));
        end
        errorbar([1,2],m,sem,'linewidth',5,'color','k');
        p = signrank(S,R);
        set(gca,'xticklabel',{'Stable','Random'});
        set(gca,'xtick',[1:2],'tickdir','out');
        xlim([0.5,2.5]);
        ylabel('Mean Interneuronal Distance [microns]');
        title(['P = ',num2str(p)]);
        
%%    
function [s,r,S,R,c] = UnpackDs(D,stabilityStatus)
    nAnimals = length(D); 
    colors = parula(nAnimals);
    
    [s,r,S,R,c] = deal([]);
    for a=1:nAnimals
        nSessions = length(D{a});
        
        for sesh=1:nSessions
            %Stable/unstable neuron indices. 
            stable = stabilityStatus.stable{a}{sesh};
            unstable = stabilityStatus.unstable{a}{sesh};
            
            %Get the mean. 
            m = nanmean(D{a}{sesh}(:));
            sd = nanstd(D{a}{sesh}(:));
            
            %Get rid of nans. 
            temp = D{a}{sesh}(stable,stable);
            temp = temp(:);
            temp(isnan(temp)) = [];
            
            %Z-score.
            %stableTemp = bsxfun(@rdivide, bsxfun(@minus, temp, m), sd);
            stableTemp = temp; 
            
            %rsamp = D{a}{sesh}(unstable,unstable);
            rsamp = randsample(union(stable,unstable),1000,true);
            rsamp = D{a}{sesh}(stable,rsamp); 
            rsamp = rsamp(:);
            rsamp(isnan(rsamp)) = [];
            
            %rsamp = bsxfun(@rdivide, bsxfun(@minus, rsamp, m), sd);
            
            
           
            s = [s; stableTemp];
            r = [r; rsamp];
            
            S = [S; mean(stableTemp)];
            R = [R; mean(rsamp)];
            
            c = [c; colors(a,:)];
        end
    end
    
end