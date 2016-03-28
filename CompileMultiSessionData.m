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
    datatypes = {   'timecells',...
                    'ratebylap',...
                    'curves',...
                    'delays',...
                    'complete',...
                    'placefields',...
                    'placefieldsnonan',...
                    'placefieldsunsmoothed',...
                    'occmaps',...
                    'placefieldpvals'};

    %Check that the arguments match template. 
    for i=1:nArgs
        assert(any(strcmp(datatypes,args{i})),['Argument ', num2str(i), ' invalid!']); 
    end
    
%% Compile. 
    for i=1:nSessions
        cd(paths{i}); 
        
        %TIME CELLS. 
        if any(strcmp('timecells',args))
            load(fullfile(pwd,'TimeCells.mat'),'TimeCells');
            DATA.timecells{i} = TimeCells; 
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
            load(fullfile(pwd,'PlaceMaps.mat'),'TMap_gauss','OccMap'); 
            for j=1:length(TMap_gauss)
                TMap_gauss{j}(OccMap==0) = 0;
            end
            DATA.placefieldsnonan{i} = TMap_gauss;  
        end
        
        %PLACE FIELDS version 3, unsmoothed. 
        if any(strcmp('placefieldsunsmoothed',args))
            load(fullfile(pwd,'PlaceMaps.mat'),'TMap_unsmoothed'); 
            DATA.placefieldsunsmoothed{i} = TMap_unsmoothed;  
        end
        
        %OCCUPANCY MAP.
        if any(strcmp('occmaps',args))
            load(fullfile(pwd,'PlaceMaps.mat'),'OccMap');
            OccMap(OccMap==0) = NaN;
            OccMap(OccMap>0) = 0;     
            DATA.occmaps{i} = OccMap;
        end
        
        %PLACE FIELD P-VALUE.
        if any(strcmp('placefieldpvals',args))
            load(fullfile(pwd,'PlaceMaps.mat'),'pval');
            DATA.placefieldpvals{i} = 1-pval;
        end
    end
    
    %Return to directory. 
    cd(initDir);
end
