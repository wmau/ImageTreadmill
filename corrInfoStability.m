function CORRS = corrInfoStability(mds,stabilityType,statType)
%corrInfoStability(mds,stabilityType,statType)
%
%   Correlate the correlation coefficient of a cell and its temporal or
%   spatial information. Takes same inputs as PartitionStats. Work on
%   commenting this.

%%
    animals = unique({mds.Animal});
    nAnimals = length(animals);
    statType = lower(statType);
    stabilityType = lower(stabilityType); 
    
%% Determine what type of neuron to analyze.
    switch stabilityType
        case 'time'
            switch statType
                case 'ti',cellGet = 'timecells'; 
                case 'si',cellGet = 'timecells'; 
                case {'fr','fluor'}, cellGet = 'timecells';
            end
        case 'place'
            switch statType
                case 'ti',cellGet = 'placecells'; 
                case 'si',cellGet = 'placecells'; 
                case {'fr','fluor'}, cellGet = 'placecells';
            end
    end
   
%% 
    [stats,corrs,neurons,lbl,preXCs,code,nonCode,d] = deal(cell(nAnimals,1)); 
    [STATS,CORRS,LABEL] = deal([]);
    for a=1:nAnimals
        %Get all the sessions for this animal.
        ssns = find(strcmp(animals{a},{mds.Animal}));
        
        [stats{a},corrs{a},neurons{a},lbl{a},preXCs{a},code{a},nonCode{a},d{a}] = deal(cell(length(ssns)-1,1));
        for s=1:length(ssns)-1
            code{a}{s} = AcquireTimePlaceCells(mds(ssns(s)),cellGet);   
            n = length(code{a}{s});
            
            switch cellGet
                case 'timecells'
                    preXCs{a}{s} = getNewTimeCells(mds(ssns(s)),mds(ssns(s+1)));
                case 'placecells'
                    preXCs{a}{s} = getNewPlaceCells(mds(ssns(s)),mds(ssns(s+1)));
            end
            
            cd(mds(ssns(s)).Location);
            load('FinalOutput.mat','NumNeurons');
            neurons{a}{s} = 1:NumNeurons;
            nonCode{a}{s} = setdiff(neurons{a}{s},union(preXCs{a}{s},code{a}{s}));
            neurons{a}{s} = EliminateUncertainMatches([mds(ssns(s)) mds(ssns(s+1))],neurons{a}{s});
                             
            cd(mds(ssns(s)).Location);
            %% Get information scores and neurons. 
            switch statType
                case 'ti'
                    load('TemporalInfo.mat','MI');
                    stats{a}{s} = MI; 
                    
                case 'si'
                    load('SpatialInfo.mat','MI');
                    stats{a}{s} = MI;       
                    
                case 'fr'
                    load('Pos_align.mat','PSAbool');
                    [n,f] = size(PSAbool);
                    d = diff([zeros(n,1) PSAbool],1,2);
                    d(d<0) = 0;
                    stats{a}{s} = sum(d,2)./f;
                    
            end  
            
            %% Correlate tuning curves. 
            switch stabilityType
                case 'time'
                    %Correlate tuning curves. 
                    corrs{a}{s} = CorrTrdmllTrace(mds(ssns(s)),mds(ssns(s+1)),neurons{a}{s});
%                     [tuningStatus,d{a}{s}] = TCRemap(mds(ssns(s)),mds(ssns(s+1)));
                    lbl{a}{s} = single(corrs{a}{s}(:,2) < .01/n); 
                    %lbl{a}{s} = single(tuningStatus(:,2) == 1); 
                    lbl{a}{s}(preXCs{a}{s}) = 2; 
                    lbl{a}{s}(nonCode{a}{s}) = 3; 
                case 'place'
                    corrs{a}{s} = CorrPlaceFields(mds(ssns(s)),mds(ssns(s+1)),neurons{a}{s});
                    %[tuningStatus,d{a}{s}] = PCRemap(mds(ssns(s)),mds(ssns(s+1)));
                    lbl{a}{s} = single(corrs{a}{s}(:,2) < .01/n); 
                    %lbl{a}{s} = single(tuningStatus(:,2) == 1); 
                    lbl{a}{s}(preXCs{a}{s}) = 2; 
                    lbl{a}{s}(nonCode{a}{s}) = 3; 
            end
            
            %% Normalize.
%             if strcmp(statType,'fr')
%                 stats{a}{s} = (stats{a}{s}-min(stats{a}{s}))./range(stats{a}{s}); 
%             else
                  stats{a}{s}(neurons{a}{s}) = zscore(stats{a}{s}(neurons{a}{s}));
%                  corrs{a}{s}(neurons{a}{s}) = zscore(corrs{a}{s}(neurons{a}{s}));
%            end

        end
        
        STATS = [STATS; cell2mat(cellfun(@(x,y) x(y),stats{a},neurons{a},'unif',0))];
        LABEL = [LABEL; cell2mat(cellfun(@(x,y) x(y),lbl{a},neurons{a},'unif',0))];
        %CORRS = [CORRS; cell2mat(cellfun(@(x,y) x(y),d{a},neurons{a},'unif',0))];
        CORRS = [CORRS; cell2mat(cellfun(@(x,y) x(y),corrs{a},neurons{a},'unif',0))];
    end
    
    CORRS(isnan(CORRS)) = 0;
    
    switch cellGet
        case 'timecells', newColor = [0 .7 .7]; 
        case 'placecells', newColor = [.78 .44 1];
    end
    stableColor = 'g';
    unstableColor = 'r';
    nonCodeColor = [.7 .7 .7];
    
    %coding = ismember(LABEL,[0,1,2]);
    coding = ismember(LABEL,[0,1]);
    [R,p] = corr(STATS(coding),CORRS(coding),'rows','complete','type','pearson');
    [sr,sp] = corr(STATS(coding),CORRS(coding),'rows','complete','type','spearman');
    [kr,kp] = corr(STATS(coding),CORRS(coding),'rows','complete','type','kendall');

    figure('Name',[statType, ' vs ',stabilityType, ' stability']);
    hold on;
    s(1) = scatter(STATS(LABEL==3),CORRS(LABEL==3),20,nonCodeColor,'filled');
    s(2) = scatter(STATS(LABEL==0),CORRS(LABEL==0),20,unstableColor,'filled');
    s(3) = scatter(STATS(LABEL==1),CORRS(LABEL==1),20,stableColor,'filled');
    s(4) = scatter(STATS(LABEL==2),CORRS(LABEL==2),20,newColor,'filled');
    yLims = get(gca,'ylim'); 
    line(median(STATS(LABEL==3))*ones(1,2),yLims,'color',nonCodeColor,'linewidth',4,'linestyle',':');
    line(median(STATS(LABEL==0))*ones(1,2),yLims,'color',unstableColor,'linewidth',4,'linestyle',':'); 
    line(median(STATS(LABEL==1))*ones(1,2),yLims,'color',stableColor,'linewidth',4,'linestyle',':'); 
    line(median(STATS(LABEL==2))*ones(1,2),yLims,'color',newColor,'linewidth',4,'linestyle',':');
    alpha(s,.4);
    
    %ylabel('Correlation Rho');
    ylabel('Correlation Coeff.'); 
    
    switch statType
        case 'ti', infoType = 'Temporal ';
        case 'si', infoType = 'Spatial ';
    end
    xlabel([infoType,'Information [bits]']);
    title({['Pearson R = ',num2str(R),', p = ',num2str(p)],...
           ['Spearman R = ',num2str(sr),', p = ',num2str(sp)],...
           ['Kendall R = ',num2str(kr),', p = ',num2str(kp)]});
    set(gca,'tickdir','out','linewidth',4,'fontsize',15);
    ylim(yLims);
    
%%
    %Check if upper and lower quartiles of correlation coefficients align
    %with information scores. 
    stable = LABEL==1;
    upperQuartile = quantile(CORRS(stable),.75); 
    lowerQuartile = quantile(CORRS(stable),.25); 
    u = CORRS > upperQuartile & stable;
    l = CORRS < lowerQuartile & stable; 
    middle = CORRS <= upperQuartile & CORRS >= lowerQuartile & stable; 
    
    [R,p] = corr(STATS(u),CORRS(u),'rows','complete','type','pearson');
    
    figure;
    scatter(STATS(u),CORRS(u),[],newColor,'filled');
    set(gca,'tickdir','out','fontsize',15);
    title({'Upper quartile',...
           ['Pearson R^{2} = ',num2str(R^2),', p = ',num2str(p)]});    
       xlabel('Information [z-bits]'); 
    ylabel('Correlation Coefficient'); 
    
    [R,p] = corr(STATS(l),CORRS(l),'rows','complete','type','pearson');
    
    figure;
    scatter(STATS(l),CORRS(l),[],newColor,'filled');
    set(gca,'tickdir','out','fontsize',15);
    title({'Lower quartile',...
           ['Pearson R^{2} = ',num2str(R^2),', p = ',num2str(p)]});        
    xlabel('Information [z-bits]'); 
    ylabel('Correlation Coefficient'); 
    
    [R,p] = corr(STATS(middle),CORRS(middle),'rows','complete','type','pearson');
    
    figure;
    scatter(STATS(middle),CORRS(middle),[],newColor,'filled');
    set(gca,'tickdir','out','fontsize',15);
    title({'Middle 50%',...
           ['Pearson R^{2} = ',num2str(R^2),', p = ',num2str(p)]});        
    xlabel('Information [z-bits]'); 
    ylabel('Correlation Coefficient'); 
    
 end