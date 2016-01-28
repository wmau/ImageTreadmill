function remap = TimeCellRemapRate(batch_session_map,base,comp,Ts)
%
%
%

%% 
    %Partition the session data. 
    MD = [base, comp];
    nSessions = length(MD); 
    dates = {MD.Date};
    sessions = [MD.Session];
      
    [TIMECELLS,~,CURVES,~,~] = CompileTimeCellData(MD,Ts);
    
%% 
    regDates = {batch_session_map.session.date};
    regSessions = [batch_session_map.session.session]; 
    
    %Eliminate first column to match indices.
    MAP = batch_session_map.map(:,2:end);  
    
    MAPcols = zeros(nSessions,1); 
    for i=1:nSessions
        try
            MAPcols(i) = find(ismember(regDates,dates{i}) & ismember(regSessions,sessions(i)));
        catch
            error(['Error in above. Possible reason: MD input has not been registered yet '...
                'for ',dates{i},' session ',num2str(sessions(i)), '. Run neuron_reg_batch...']);
        end
       
    end
   
    %Rows on MAP corresponding to where time cells are found in the base
    %session. 
    MAProws = find(ismember(MAP(:,MAPcols(1)),TIMECELLS{1})); 
    
%%
    %Find time resolution of the tuning curves. 
    tResolution = zeros(nSessions,1);
    for thisSession=1:nSessions
        nBins = length(CURVES{thisSession}.tuning{1});
        tResolution(thisSession) = Ts(thisSession)/nBins;
    end
    
    %If there are multiple time resolutions, throw an error.
    assert(length(unique(tResolution))==1,'Different time resolutions found in tuning curves! Rerun FindTimeCells!');
    tResolution = unique(tResolution);      %seconds. 
    
    window = 0.5;                                   %seconds.
    binWindow = round(1/tResolution*window);        %bins.
    remap = nan(max(MAProws),nSessions); 
    
    for thisRow=MAProws'
        neurons = MAP(thisRow,MAPcols);
        baseSig = CURVES{1}.sig{neurons(1)}; 
        
        %Find humps in the significance curve (essentially translates to
        %the tuning curve).
        ccBase = bwconncomp(baseSig);
        nHumps = length(ccBase.PixelIdxList); 
        
        %Get bin numbers corresponding to the middle of the hump. .
        baseBins = zeros(nHumps,1); 
        for i=1:nHumps
            baseBins(i) = round(mean(ccBase.PixelIdxList{i}));
        end
        
        for thisSession=2:nSessions
            if ismember(neurons(thisSession),TIMECELLS{thisSession})
                compSig = CURVES{thisSession}.sig{neurons(thisSession)};
                
                %Find humps in the significance curve of the comparison
                %session. 
                ccComp = bwconncomp(compSig);
                nHumps = length(ccComp.PixelIdxList); 
                
                %Get bin numbers corresponding to the middle of the hump. 
                compBins = zeros(nHumps,1);
                for i=1:nHumps
                    compBins(i) = round(mean(ccComp.PixelIdxList{i}));
                end
                
                %All possible combinations. 
                [b,c] = meshgrid(baseBins,compBins); 
                
                %Subtract.
                if any(abs(b-c) <= binWindow)
                    remap(thisRow,thisSession) = 0; 
                else
                    remap(thisRow,thisSession) = 1; 
                end
                
            end
            
        end
        
    end
    
    keyboard;
    
end