function DATA = CompileMultiSessionData(MD,args)
%DATA = CompileMultiSessionData(MD,args)
%
%   Compile data from multiple sessions. Data collected is dependent on the
%   content of args. 
%
%   INPUTS
%       MD: Entries indicating what session. 
%
%       args: Cell array containing strings that could be either...
%           timecells - time cell indices.
%           ratebylap - rasters.
%           curves - tuning, shuffled curves along with significance and
%               p-values, confidence intervals. 
%           delays - delay durations.
%           complete - whether the run was complete or not 
%           placefields - place fields. 
%           occmaps - occupancy maps. 
%           placefieldpvals - place field entropy p-values. 
%
%   OUTPUT
%       DATA: structure array containing fields for each arg. Each field is
%           a Sx1 cell array for each session. 
%

%% Preliminary steps. 
    initDir = pwd; 
    
    %Partition the session data.
    nSessions = length(MD);
    paths = {MD.Location}; 

    nArgs = length(args); 
    args = lower(args); 
    
    %Template. 
    datatypes = {   'ft',...
                    'ttl',...
                    't',...
                    'timecells',...
                    'ratebylap',...
                    'curves',...
                    'delays',...
                    'complete',...
                    'placefields',...
                    'placefieldsnonan',...
                    'placefieldsunsmoothed',...
                    'runoccmaps',...
                    'placefieldpvals',...
                    'a',...
                    'rawtrdmll',...
                    'difftrdmll',...
                    'tracetrdmll',...
                    'si',...
                    'ti',...
                    'placecells'};

    %Check that the arguments match template. 
    for i=1:nArgs
        assert(any(strcmp(datatypes,args{i})),['Argument ', num2str(i), ' invalid!']); 
    end
    
%% Compile. 
    for i=1:nSessions
        cd(paths{i}); 
        %FT.
        if any(strcmp('ft',args))
            load(fullfile(pwd,'Pos_align.mat'),'FT');
            DATA.ft{i} = FT;
        end
        
        %TODAYTREADMILLLOG.
        if any(strcmp('ttl',args))
            load(fullfile(pwd,'TimeCells.mat'),'TodayTreadmillLog');
            DATA.ttl{i} = TodayTreadmillLog;
        end
        
        %T.
        if any(strcmp('t',args))
            load(fullfile(pwd,'TimeCells.mat'),'T');
            DATA.t{i} = T;
        end
        
        %TIME CELLS. 
        if any(strcmp('timecells',args))
            load(fullfile(pwd,'TimeCells.mat'),'TimeCells');
            load(fullfile(pwd,'TemporalInfo.mat'),'sig');
            DATA.timecells{i} = intersect(find(sig),TimeCells); 
            %DATA.timecells{i} = TimeCells;
        end
        
        %RASTERS.
        if any(strcmp('ratebylap',args))
            load(fullfile(pwd,'TimeCells.mat'),'ratebylap'); 
            DATA.ratebylap{i} = ratebylap; 
        end
        
        %CURVES.
        if any(strcmp('curves',args))
            load(fullfile(pwd,'TimeCells.mat'),'curves'); 
            DATA.curves{i} = curves; 
        end
        
        %DELAYS.
        if any(strcmp('delays',args))
            load(fullfile(pwd,'TimeCells.mat'),'TodayTreadmillLog'); 
            DATA.delays{i} = TodayTreadmillLog.delaysetting; 
        end
        
        %COMPLETE. 
        if any(strcmp('complete',args))
            load(fullfile(pwd,'TimeCells.mat'),'TodayTreadmillLog'); 
            DATA.complete{i} = logical(TodayTreadmillLog.complete);  
        end
        
        %PLACE FIELDS. 
        if any(strcmp('placefields',args))
            load(fullfile(pwd,'PlaceMaps.mat'),'TMap_gauss','OccMap'); 
            for j=1:length(TMap_gauss)
                TMap_gauss{j}(OccMap==0) = nan;
            end
            DATA.placefields{i} = TMap_gauss;  
        end
        
        %PLACE FIELDS version 2, no NaNs. 
        if any(strcmp('placefieldsnonan',args))
            load(fullfile(pwd,'PlaceMaps.mat'),'TMap_gauss','RunOccMap'); 
            for j=1:length(TMap_gauss)
                TMap_gauss{j}(RunOccMap==0) = 0;
            end
            DATA.placefieldsnonan{i} = TMap_gauss;  
        end
        
        %PLACE FIELDS version 3, unsmoothed. 
        if any(strcmp('placefieldsunsmoothed',args))
            load(fullfile(pwd,'PlaceMaps.mat'),'TMap_unsmoothed'); 
            DATA.placefieldsunsmoothed{i} = TMap_unsmoothed;  
        end
        
        %OCCUPANCY MAP.
        if any(strcmp('runoccmaps',args))
            load(fullfile(pwd,'PlaceMaps.mat'),'RunOccMap');
            RunOccMap(RunOccMap==0) = NaN;
            RunOccMap(RunOccMap>0) = 0;     
            DATA.runoccmaps{i} = RunOccMap;
        end
        
        %PLACE FIELD P-VALUE.
        if any(strcmp('placefieldpvals',args))
            load(fullfile(pwd,'Placefields.mat'),'pval');
            DATA.placefieldpvals{i} = pval;
        end
        
        %ADJACENCY MATRIX.
        if any(strcmp('a',args))
            load(fullfile(pwd,'graphData_p.mat'),'A');
            DATA.A{i} = A; 
        end
        
        %RAW TRACES ON TREADMILL.
        if any(strcmp('rawtrdmll',args))
            load(fullfile(pwd,'TreadmillTraces.mat'),'rawtrdmll');
            DATA.rawtrdmll{i} = rawtrdmll; 
        end
        
        %DIFFERENTIAL TRACES ON TREADMILL.
        if any(strcmp('difftrdmll',args))
            load(fullfile(pwd,'TreadmillTraces.mat'),'difftrdmll');
            DATA.difftrdmll{i} = difftrdmll; 
        end
        
        %Z-SCORE TRACES ON TREADMILL.
        if any(strcmp('tracetrdmll',args))
            load(fullfile(pwd,'TreadmillTraces.mat'),'tracetrdmll');
            DATA.tracetrdmll{i} = tracetrdmll; 
        end
        
        %SPATIAL INFORMATION
        if any(strcmp('si',args))
            load(fullfile(pwd,'SpatialInfo.mat'),'MI');
            DATA.si{i} = MI; 
        end
        
        %TEMPORAL INFORMATION.
        if any(strcmp('ti',args))
            load(fullfile(pwd,'TemporalInfo.mat'),'MI');
            DATA.ti{i} = MI; 
        end
        
        %PLACE CELLS.
        if any(strcmp('placecells',args))
            load(fullfile(pwd,'Placefields.mat'),'pval');
            load(fullfile(pwd,'PlacefieldStats.mat'),'PFnHits','bestPF');
            idx = sub2ind(size(PFnHits), 1:size(PFnHits,1), bestPF');
            DATA.placecells{i} = find(pval < .01 & PFnHits(idx) > 4); 
        end
    end
    
    %Return to directory. 
    cd(initDir);
end
