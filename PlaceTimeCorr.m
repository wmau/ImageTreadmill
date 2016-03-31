function [pCorr,tCorr,MAP,MAPcols,DATA,noi] = PlaceTimeCorr(MAPMD,MD1,MD2,noi)
%[pCorr,tCorr,MAP,MAPcols,DATA,noi] = PlaceTimeCorr(MAPMD,MD1,MD2)
%
%   Computes the correlation of the place fields of time cells as well as
%   the correlation of the temporal tuning curve (all unsmoothed). 
%
%   INPUTS
%       mapMD: MD entry of where batch_session_map lives. 
%
%       MD1: Base MD entry. The function will look at time cells here. 
%
%       MD2: Comparison MD entry. Correlations will be done to this
%       session.
%
%   OUTPUTS
%       pCorr: Nx2 matrix of correlation results. First column is rho,
%       second column is p-value. 
%
%       tCorr: Same as pCorr, but for temporal tuning curve. 
%
%       MAP: batch_session_map.
%
%       MAPcols: Columns corresponding to the dates specified. 
%
%       DATA: Data compiled over multiple sessions. See
%       CompileMultiSessionData. 
%
%       noi: Neurons of interest, indexes MD1.  
%

%% Compile data across multiple sessions. 
    sessions = [MD1,MD2]; 
    dates = {sessions.Date};
    sessionNums = [sessions.Session];
    
    DATA = CompileMultiSessionData(sessions,...
        {'timecells','curves','placefields','placefieldsunsmoothed',...
        'placefieldpvals','ratebylap','delays','complete','occmaps'}); 
   
    %Neurons of interest, aka time cells with place fields.
    nNeurons = length(DATA.placefieldsunsmoothed{1}); 
    
    %Load neuron mapping.
    load(fullfile(MAPMD.Location,'batch_session_map.mat'));
    [MAP,MAPcols] = FilterMAPDates(batch_session_map,dates,sessionNums);
        
    MAProws = find(ismember(MAP(:,MAPcols(1)),noi));        %Indices for neurons of interest. 
    goodrows = MAProws(~any(MAP(MAProws,MAPcols)==0,2));    %Cut out neurons that weren't mapped in session 2. 
    noi = noi(ismember(noi,MAP(goodrows,MAPcols(1))));      %Exclude unmapped neurons. .

    %Preallocate. 
    pCorr = nan(nNeurons,2); 
    tCorr = nan(nNeurons,2);    
    for n1=1:nNeurons
        if ismember(n1,noi)
            n2 = MAP(MAP(:,MAPcols(1))==n1,MAPcols(2));
            
            %Place fields.
            pf1 = DATA.placefieldsunsmoothed{1}{n1}(:);
            pf2 = DATA.placefieldsunsmoothed{2}{n2}(:);          
            [pCorr(n1,1),pCorr(n1,2)] = corr(pf1,pf2);
            
            %Time fields.
            tf1 = DATA.curves{1}.tuning{n1}; 
            tf2 = DATA.curves{2}.tuning{n2};   
            
            [tCorr(n1,1),tCorr(n1,2)] = corr(tf1',tf2');             
        end
    end
    
end