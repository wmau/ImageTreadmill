function xCorrTrdmllTraces(md,tracetype,corrType)
%
%
%

%% Set up.
    tic; 
    cd(md.Location);
    load('Pos_align.mat',tracetype);
    switch tracetype
        case 'rawtrace'
            FT = rawtrace; clear rawtrace; 
        case 'difftrace' 
            FT = difftrace; clear difftrace;
    end
      
    %Normalize by max.
    m = max(FT,[],2);
    FT = FT./repmat(m,1,size(FT,2));
    
    %Trim treadmill indices.
    load('TimeCells.mat','T','TodayTreadmillLog'); 
    inds = TrimTrdmllInds(TodayTreadmillLog,T);
    
    %Get active neurons.
    nNeurons = size(FT,1);
    [~,active] = nNeuronsActiveonTM(md);
    
%% Build rasters.
    rasters = cell(1,nNeurons);
    for n=active
        rasters{n} = buildRasterTrace(inds,FT,n);
    end

    nComparisons = nNeurons*nNeurons;
%% Perform cross-correlations and permutation tests. 
    B = 500;
    R = cell(nNeurons);
     
    nLaps = size(rasters{active(1)},1); 
    resolution = 2;
    updateInc = round(nComparisons/(100/resolution));
    p = ProgressBar(100/resolution);
    parpool('local');
    parfor c=1:nComparisons
        [src,snk] = ind2sub([nNeurons,nNeurons],c); 
        
%% Time shuffle. 
        if src<snk && ismember(src,active) && ismember(snk,active)
            [R{c}.trials,lags] = xcorr_by_laps(rasters{src},rasters{snk},corrType);
            R{c}.curve = mean(R{c}.trials);
            R{c}.CI = 1.96*std(R{c}.trials)/sqrt(nLaps);
            
            lLags = length(lags);
            tshffcurves = zeros(B,lLags);
            for i=1:B
                %Shuffle time.
                tempsrc = rasters{src};
                tempsrc = permuteTime(tempsrc); 
                
                %Cross correlation then take the mean. 
                tshfftrials = xcorr_by_laps(tempsrc,rasters{snk},corrType);
                tshffcurves(i,:) = mean(tshfftrials);
            end
            
            %Sort for CIs.
            tshffcurves = sort(tshffcurves); 
            
            %Means of all curves and their confidence intervals.
            R{c}.tshff.mu = mean(tshffcurves);
            R{c}.tshff.upper = tshffcurves(round(.95*B),:);
            R{c}.tshff.lower = tshffcurves(round(.05*B),:); 
            
            %Significance. 
            R{c}.sigt = sum(repmat(R{c}.curve,B,1) <= tshffcurves)./B;
            
%% If significant difference found by shuffling time...
            if any(R{c}.sigt < .05)
                %Shuffle trial identity.
                trlshffcurves = zeros(B,lLags);
                
                for i=1:B
                    tempsrc = rasters{src}(randperm(nLaps),:);

                    %Cross correlation.
                    trlshfftrials = xcorr_by_laps(tempsrc,rasters{snk},corrType);
                    trlshffcurves(i,:) = mean(trlshfftrials);
                end
                
                %Sort for CIs.
                trlshffcurves = sort(trlshffcurves); 

                %Means of all curves and their confidence intervals.
                R{c}.trlshff.mu = mean(trlshffcurves);
                R{c}.trlshff.upper = trlshffcurves(round(.95*B),:);
                R{c}.trlshff.lower = trlshffcurves(round(.05*B),:); 

                %Significance. 
                R{c}.sigtrls = sum(repmat(R{c}.curve,B,1) <= trlshffcurves)./B;
                
                %Final significance logical.
                R{c}.sig = R{c}.sigtrls < .05 & R{c}.sigt < .05 & ...
                    (R{c}.curve - R{c}.CI) > R{c}.trlshff.upper & ...
                    (R{c}.curve - R{c}.CI) > R{c}.tshff.upper;
            else 
                R{c}.sig = false(size(R{c}.sigt));
                R{c}.sigtrls = false(size(R{c}.sigt));
            end
        end
        
        if round(c/updateInc) == (c/updateInc)
            p.progress;     
        end
    end
    p.stop;
    delete(gcp); 
    
%% Build adjacency matrix. 
    %Get the lag vector again. 
    [~,lags] = xcorr_by_laps(rasters{active(1)},rasters{active(1)},corrType);
    lLags = length(lags);
    
    %Time vector and preallocate. 
    t = linspace(-T,T,lLags); 
    A = zeros(nNeurons);
    window = 2; 
    for i=active   
        nonI = active(active>i);
        for j=nonI
            sigs = R{i,j}.sig;      %Package.
            %tuningcurve = R.smoothed{i,j};  %Get the tuning curve so you can look for peak.
            %tuningcurve(~sigs) = nan;       %Only look at significant points. 

            %Define time window. 
            %tuningcurve(t < -window | t > window) = nan;

            %[~,peak] = max(tuningcurve); 
            peak = mean(t(sigs));

            %Place edges. 
            if peak < 0 && peak > -window
                A(i,j) = 1;                 %If peak xcorr is below 0, i connects to j.
            elseif peak > 0 && peak < window
                A(j,i) = 1;                 %If peak xcorr is above 0, j connects to i. 
            elseif peak == 0
                A(i,j) = 1;
                A(j,i) = 1; 
            end

        end
    end
    
    %Divide by frame rate.
    lags = lags./20;
    
%% Excise overlap.
    [edgesRemoved,A] = ExciseOverlap(md,A);
  
%% Save.
    elapsed = toc;   
    save('XCorr.mat','R','A','lags','edgesRemoved','elapsed',...
        'corrType','tracetype','-v7.3');
    
end