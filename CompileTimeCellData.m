function [TIMECELLS,RATEBYLAP,CURVES,DELAYS,COMPLETE] = CompileTimeCellData(MD,Ts)
%[TIMECELLS,RATEBYLAP,CURVES,DELAYS,COMPLETE] = CompileTimeCellData(MD)
%
%   Extracts data from saved TimeCells.mat for all the entries listed in
%   the input variable MD. If not found, will run FindTimeCells()
%   automatically.
%
%   INPUT
%       MD: MasterDirectory entry. 
%
%       Ts: Vector, same length as MD specifying delay durations you want
%       to look at. 
%
%   OUTPUTS
%   All outputs are cell arrays the same size as MD. Below describes the
%   contents found within a single cell referring to a particular session.
%
%       TIMECELLS: Vector containing indices referencing FT.
%
%       RATEBYLAP: LxTxN matrix (L=laps, T=time bins, N=neurons). These are
%       complete ratebylap matrices, containing incomplete laps and laps
%       that were running at different durations. 
%       
%       CURVES: Structure array with fields...
%       Each field contains a Nx1 cell. 
%           tuning: 1xT vector, trial average of ratebylap for that neuron.
%           
%           shuffle: BxT matrix (B=shuffles), permuted responses. 
%
%           sig: 1xT logical, significant or not. 
%
%           p: 1xT vector, p values. 
%
%           ci: 2xT vector, 95% confidence intervals. First row is upper
%           bound. 
%
%       DELAYS: Lx1 vector containing delay durations. 
%
%       COMPLETE: Lx1 logical, completed lap or not (fulfilling all
%       criteria such as whether mouse was on treadmill and whether the run
%       was complete). 
%   

%% Organization. 
    assert(length(Ts)==length(MD),'Length of MD and Ts must be the same!'); 

    initDir = pwd; 

    %Partition the session data.
    nSessions = length(MD);
    dates = {MD.Date}; 
    paths = {MD.Location}; 
    animals = {MD.Animal};
    sessions = [MD.Session];
    
   %Preallocate.
    TIMECELLS = cell(nSessions,1);
    RATEBYLAP = cell(nSessions,1);
    CURVES = cell(nSessions,1);
    DELAYS = cell(nSessions,1); 
    COMPLETE = cell(nSessions,1); 
  
%% Gather all time cell data. 
    for i=1:nSessions
        cd(paths{i});
        
        try 
            load(fullfile(paths{i},'TimeCells.mat'),'TimeCells','ratebylap','TodayTreadmillLog','curves','T');
            
            %Throw an error if you specify a T that doesn't match the one
            %that you already ran. 
            if T~=Ts(i)
                disp(['Delay duration specified in T(',num2str(i),') is different '...
                    'from the one saved for ',dates{i},'!']);
                
                return; 
            end
            
            %Archive. 
            TIMECELLS{i} = TimeCells;
            RATEBYLAP{i} = ratebylap; 
            CURVES{i} = curves; 
            DELAYS{i} = TodayTreadmillLog.delaysetting; 
            COMPLETE{i} = logical(TodayTreadmillLog.complete); 
            
        catch
            [TIMECELLS{i},RATEBYLAP{i},CURVES{i},~,~,TodayTreadmillLog]...
                = FindTimeCells(animals{i},dates{i},sessions(i),Ts(i));
            
            DELAYS{i} = TodayTreadmillLog.delaysetting; 
            COMPLETE{i} = logical(TodayTreadmillLog.complete);
        end
        
    end
    
    cd(initDir);
end