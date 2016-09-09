function PlotStability(mds)
%
%
%

%% 
    %What animals are in md?
    animals = unique({mds.Animal}); 
    nAnimals = length(animals); 
    
    DATA = CompileMultiSessionData(mds,{'t'});
    
    minDaysFromRef = 0;
    maxDaysFromRef = 0;
    %Get mapMD for each animal.
    for a=1:nAnimals
        %Get all the MDs for that animal and their corresponding mapMDs. 
        mdIdx = ismember({mds.Animal},animals{a});
        theseMDs = mds(mdIdx);
        mapMDs(a) = getMapMD(theseMDs);        
        
        %Get the middle MD. 
        middle = ceil(length(theseMDs)/2);
        baseMD = theseMDs(middle);
        
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
    stability =     cell(1,nAnimals); 
    meanStability = cell(1,nAnimals); 
    dayrank =       cell(1,nAnimals); 
    sStability =    cell(1,nAnimals); 
    STABILITY = 	cell(1,maxnSessions); 
    sSTABILITY =    cell(1,maxnSessions); 
    vStability = [];
    vDays = [];
    B = 1000;
    
    for a=1:nAnimals
        %Get all the MDs for that animal.
        mdIdx = ismember({mds.Animal},animals{a});
        theseMDs = mds(mdIdx);
        
        %Make sure the dates are sorted chronologically. 
        dates = {theseMDs.Date};
        theseMDs = sortbyDate(dates,theseMDs); 
        
        %Get the middle MDs and the rest. 
        middle = ceil(length(theseMDs)/2);
        baseMD = theseMDs(middle);
        compMDs = theseMDs([1:middle-1,middle+1:end]);
        
        %Get neurons across days and rank TCs according to the base
        %session.
        Ts = [DATA.t{mdIdx}];
        [~,sortedPeaks] = multiPastalkovaPlot(mapMDs(a),baseMD,compMDs,Ts,false);
        
        %Compute stability. 
        stability{a} = sDistanceMetric(sortedPeaks);
        
        %Shuffle the peaks before sorting by date. 
        shuffledPeaks = nan(size(sortedPeaks));             %Preallocate for shuffle test.
        shuffledPeaks(:,1) = sortedPeaks(:,1);              %Don't shuffle base session. 
        
        %Sort by date then append onto long vector. 
        stability{a} = sortbyDate({baseMD.Date compMDs.Date},stability{a});
        
        %Do same test on shuffled data. 
        [nTCs,nSessions] = size(stability{a});
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
        
        %Take the mean of the control and the empirical distribution. 
        sMeans = sort(squeeze(nanmean(sStability{a}))');
        meanStability{a} = nanmean(stability{a}); 
        
        dayrank{a} = zeros(1,nSessions);     
        for s=1:length(theseMDs)
            %Get days apart.  
            reltoRef = daysact(baseMD.Date,theseMDs(s).Date); 
            dayrank{a}(s) = reltoRef; 
            [~,temp] = ismember(reltoRef,DaysFromRef); 
            STABILITY{temp} = [STABILITY{temp}; stability{a}(:,s)];
            sSTABILITY{temp} = [sSTABILITY{temp}; sMeans(:,s)];
            
            vStability = [vStability; stability{a}(:,s)];
            vDays = [vDays; reltoRef*ones(nTCs,1)];
        end
        
        %Plot animal points. 
        scatter(dayrank{a},meanStability{a},'d'); hold on; 
    end
    
    %Calculate shuffle control mean and CI.
    ci = nan(maxnSessions,2);
    for s=1:maxnSessions
        N = length(sSTABILITY{s});
        ci(s,1) = sSTABILITY{s}(round(0.05*N));
        ci(s,2) = sSTABILITY{s}(round(0.95*N));
    end
    sMEANSTABILITY = cellfun(@nanmean,sSTABILITY);  %Shuffle mean.
    ci = abs(ci-repmat(sMEANSTABILITY',1,2));       %Shuffle CI.
    MEANSTABILITY = cellfun(@nanmean,STABILITY);    %Empirical mean.

    %Plot.
    %plot(DaysFromRef,MEANSTABILITY,'r','linewidth',2); 
    t = minDaysFromRef:0.01:maxDaysFromRef; 
    p = pchip(DaysFromRef,MEANSTABILITY,t);
    plot(t,p,'r','linewidth',2); 
    l = boundedline(DaysFromRef,sMEANSTABILITY,ci,'alpha');
        l.Color = [.5 .5 1]; l.LineStyle = '--'; 
        xlabel('Days from Reference'); 
        ylabel('Stability Score'); 
        set(gca,'ticklength',[0 0]); 
     
    keyboard;
    
end