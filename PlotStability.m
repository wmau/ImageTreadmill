function PlotStability(mds,type,mid)
%PlotStability(mds,type,mid)
%
%   Plot the stability of time cells over days for multiple animals. You
%   can use different metrics for stability. Currently, the two options are
%   'distance' and 'correlation'. 'Distance' takes the base session and
%   draws a line down the sorted peaks of time cell tuning curves then
%   examines other sessions' tuning curve peaks and their relative distance
%   to the corresponding point on the line. 'Correlation' does pairwise
%   Spearman correlations for each time cell from the base session to all
%   other sessions. 
%
%   INPUTS
%       mds: List of session entries from multiple animals. The function
%       will automatically group sessions from the same animal and sort
%       chronologically. 
%
%       type
%           'distance': Deviation of sorted peak location from original
%           base sesssion.
%           'correlation': Pairwise Spearman correlations for each time
%           cell. 
%       
%       mid: Logical. If true, makes the reference the middle session (if
%       even number of sessions, takes the earlier of the two middle
%       sessions). 
%

%% Set up.
    %What animals are in md?
    animals = unique({mds.Animal}); 
    nAnimals = length(animals); 
    colors = parula(nAnimals);
    
    %Compile data from each session. 
    DATA = CompileMultiSessionData(mds,{'t'});
    
    %We want to find the farthest backward and forward dates relative to
    %the reference session. Set these to 0 initially. 
    minDaysFromRef = 0; 
    maxDaysFromRef = 0;
    %Get folder housing neuron mappings for each animal.
    for a=1:nAnimals
        %Get all the MDs for that animal and their corresponding mapMDs. 
        mdIdx = ismember({mds.Animal},animals{a});
        theseMDs = mds(mdIdx);
        mapMDs(a) = getMapMD(theseMDs);        
        
        %If mid, get the middle session for each animal. Otherwise, first. 
        if mid, baseIdx = ceil(length(theseMDs)/2);
        else baseIdx = 1; end
        baseMD = theseMDs(baseIdx);
        
        %Find the day spread relative to the middle. 
        for s=1:length(theseMDs)
            %Calculate how many days apart this session is from the
            %reference.
            daysApart = daysact(baseMD.Date,theseMDs(s).Date);
            
            %Find the farthest day going back in time. 
            if daysApart < minDaysFromRef
                minDaysFromRef = daysApart;
            %And forward. 
            elseif daysApart > maxDaysFromRef 
                maxDaysFromRef = daysApart;
            end
        end
    end
    %Make vector from -X : +Y days from reference. 
    DaysFromRef = minDaysFromRef:maxDaysFromRef; 
    maxnSessions = length(DaysFromRef); 
    
    %Preallocate things. 
    stability =     cell(1,nAnimals);       %Matrix of stability scores per time cell per session per animal. 
    meanStability = cell(1,nAnimals);       %Vector of stability scores per animal per session.
    dayrank =       cell(1,nAnimals);       %Range of dates from -X to +Y per animal.
    sStability =    cell(1,nAnimals);       %TC x Sess x B matrix of stability scores per animal.
    STABILITY = 	cell(1,maxnSessions);   %Compiled stability scores, removing animal identity, only keeping date range. 
    sSTABILITY =    cell(1,maxnSessions);   %Same as above, but for shuffled data. 
    B = 500;
    
    figure; hold on; 
    for a=1:nAnimals
        %Get all the MDs for that animal.
        mdIdx = ismember({mds.Animal},animals{a});
        theseMDs = mds(mdIdx);
        
        %Make sure the dates are sorted chronologically. 
        dates = {theseMDs.Date};
        theseMDs = sortbyDate(dates,theseMDs); 
        
        %Get the reference session and the comparators. 
        if mid      %Middle sessions per animal.
            baseIdx = ceil(length(theseMDs)/2);
        else        %First session per animal.
            baseIdx = 1; 
        end
        baseMD = theseMDs(baseIdx);
        compMDs = theseMDs([1:baseIdx-1,baseIdx+1:end]);
        
        %Get neurons across days and rank TCs according to the base
        %session.
        Ts = [DATA.t{mdIdx}];
        [normtilemat,sortedPeaks] = msPastalkovaPlot(mapMDs(a),baseMD,compMDs,Ts,false);
        [nTCs,nSessions] = size(sortedPeaks);
        
        %Compute stability. 
        switch type
            case 'distance'     %Distance-from-line metric. 
                stability{a} = sDistanceMetric(sortedPeaks);
                
                %Shuffle the peaks before sorting by date. 
                shuffledPeaks = nan(size(sortedPeaks));             %Preallocate for shuffle test.
                shuffledPeaks(:,1) = sortedPeaks(:,1);              %Don't shuffle base session. 
                           
            case 'correlation'  %Correlation metric. 
                stability{a} = sCorrMetric(normtilemat);
                
                %Initialize.  
                shuffledMat = cell(1,nSessions); 
                shuffledMat{1} = normtilemat{1};
                
        end
        
        %Sort by date. 
        stability{a} = sortbyDate({baseMD.Date compMDs.Date},stability{a});
        
        %Generate surrogate data. 
        switch type
            case 'distance'
                sStability{a} = nan([size(sortedPeaks),B]); 
                
                for i=1:B
                    for s=2:nSessions
                        %Shuffle ranks.
                        shuffledPeaks(:,s) = sortedPeaks(randperm(nTCs),s);
                    end
                
                %Get stability metric and sort by date. 
                sStability{a}(:,:,i) = sDistanceMetric(shuffledPeaks);
                sStability{a}(:,:,i) = sortbyDate({baseMD.Date compMDs.Date},sStability{a}(:,:,i));
               
                end
            case 'correlation'
                for i=1:B
                    for s=2:nSessions
                        %Shuffle ranks. 
                        shuffledMat{s} = normtilemat{s}(randperm(nTCs),:);
                    end
                    
                    %Correlate time cells then sort by date. 
                    sStability{a}(:,:,i) = sCorrMetric(shuffledMat);
                    sStability{a}(:,:,i) = sortbyDate({baseMD.Date compMDs.Date},sStability{a}(:,:,i));
                end
        end
        
        %Take the mean of the control and the empirical distribution. 
        sMeans = sort(squeeze(nanmean(sStability{a}))');
        meanStability{a} = nanmean(stability{a}); 
        
        dayrank{a} = zeros(1,nSessions);     
        for s=1:length(theseMDs)
            %Get days apart.  
            reltoRef = daysact(baseMD.Date,theseMDs(s).Date); 
            dayrank{a}(s) = reltoRef; 
            [~,temp] = ismember(reltoRef,DaysFromRef);  %Find appropriate relative date.
            STABILITY{temp} = [STABILITY{temp}; stability{a}(:,s)]; %Append. 
            sSTABILITY{temp} = [sSTABILITY{temp}; sMeans(:,s)];        
        end
        
        %Plot animal points. 
        scat = scatter(dayrank{a},meanStability{a},80,colors(a,:),'d'); alpha(scat,0.5);
    end
    
    %Calculate shuffle control mean and CI.
    ci = nan(maxnSessions,2);
    for s=1:maxnSessions
        N = length(sSTABILITY{s});
        ci(s,1) = sSTABILITY{s}(round(0.05*N));     %Lower bound.
        ci(s,2) = sSTABILITY{s}(round(0.95*N));     %Upper bound. 
    end
    sMEANSTABILITY = cellfun(@nanmean,sSTABILITY);  %Shuffle mean.
    ci = abs(ci-repmat(sMEANSTABILITY',1,2));       %Shuffle CI.
    MEANSTABILITY = cellfun(@nanmean,STABILITY);    %Empirical mean.

    %Not smoothed. 
    %plot(DaysFromRef,MEANSTABILITY,'r','linewidth',2); 
    
    %Smooth vectors.  
    smoothDays = minDaysFromRef:0.01:maxDaysFromRef; 
    smoothMean = pchip(DaysFromRef,MEANSTABILITY,smoothDays);
    smoothShuffMean = pchip(DaysFromRef,sMEANSTABILITY,smoothDays); %Mean shuffled. 
    smoothCI(:,1) = pchip(DaysFromRef,ci(:,1),smoothDays);          %Lower bound CI.
    smoothCI(:,2) = pchip(DaysFromRef,ci(:,2),smoothDays);          %Upper bound CI.
    
    %Plot.
    plot(smoothDays,smoothMean,'k','linewidth',2); 
    omit = smoothDays > -1 & smoothDays < 1;    
    smoothShuffMean(omit) = nan; 
    smoothCI(omit,:) = nan;
    l = boundedline(smoothDays,smoothShuffMean,smoothCI,'alpha','nan','gap');
        l.Color = [.5 .5 1]; l.LineStyle = '--'; 
        xlabel('Days from Reference'); 
        ylabel('Stability Score'); 
        set(gca,'ticklength',[0 0]); 

end