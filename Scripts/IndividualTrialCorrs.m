%% Set up.
    clear;
    
    %Load MD entries. 
    loadMD;
    %MD(292:295) = G45.
    %MD(296:299) = G48.
    %MD(300:304) = Bellatrix.
    %MD(305:309) = Polaris.
    fulldataset = MD(292:309);
    
    animals = unique({fulldataset.Animal});
    nAnimals = length(animals);
    
    %Parameters to change.
    codingCells = 'timecells';      %Options: timecells or placecells
    z = true;                       
    similarityMetric = 'corr';      %Options: corr or innerproduct.
    
    switch similarityMetric 
        case 'corr'
            yLabel = 'Mean similarity (R)';
        case 'innerproduct'
            yLabel = 'Mean similarity (inner product)';
    end
    
    %Preallocate.
    ALL_Rs = [];
    for a=1:nAnimals
        %Get relevant sessions.
        ssns = find(strcmp(animals{a},{fulldataset.Animal}));
        
        %Do PV correlations.
        [R,lapNum,sessionNum] = PVTrialCorr2(fulldataset(ssns),'similarityMetric',...
            similarityMetric,'codingCells',codingCells);
        R_all = nanmean(R,3);
        [~,sessionRs] = avgCorrMatrixOverAllDays(R_all,lapNum,sessionNum,'binByNTrials');
        
        ALL_Rs = [ALL_Rs; sessionRs];
    end
    
    matSize = min(cellfun('length',ALL_Rs));
    
    nSessions = length(ALL_Rs);
    R = zeros(matSize,matSize,nSessions);
    for s=1:nSessions
        R(:,:,s) = ALL_Rs{s}(1:matSize,1:matSize);
    end
    
    figure('position',[360 170 765 570]);
    subplot(1,2,1);
    imagesc(nanmean(R,3)); 
    colormap hot; 
    axis equal; axis tight; 
    set(gca,'ticklength',[0 0],'linewidth',4,'fontsize',12); 
    xlabel('Trials'); ylabel('Trials');
    
    [m,sem,diags] = collapseByLag(R); 
    subplot(1,2,2);
    errorbar(0:matSize-1,m,sem,'linewidth',4,'capsize',0);
    set(gca,'tickdir','out','linewidth',4,'fontsize',12);
    xlabel('Lag','fontsize',15); 
    ylabel(yLabel,'fontsize',15);
    xlim([-.5 matSize+.5]);
    
    X = cell2mat(diags');
    lags = [];
    for i=1:length(diags)
        lags = [lags (i-1).*ones(1,length(diags{i}))];
    end
    [~,~,stats] = anovan(X,{lags},'display','off');