function percentStable = multiProportionStable(mapMD,MD1,MDs,neurontype)
%percentStable = multiProportionStable(mapMD,MD1,MDs,neurontype)
%
%   Simple function to quickly do correlations of place or time cells
%   across sessions. 
%
%   INPUTS
%       mapMD: MD entry where batch_session_map lives. 
%
%       MD1: 

%%
    nCompSessions = length(MDs);
    load(fullfile(MD1.Location,'ProcOut.mat'),'NumNeurons'); 
    neurontype = lower(neurontype); 
    r = zeros(NumNeurons,2,nCompSessions); 
    switch neurontype
        case 'time'
            load(fullfile(MD1.Location,'TimeCells.mat'),'TimeCells'); 
            n = length(TimeCells);
            
            for s=1:nCompSessions
                [~,r(:,:,s)] = PlaceTimeCorr(mapMD,MD1,MDs(s),TimeCells); 
            end
        case 'place'
            load(fullfile(MD1.Location,'PlaceMaps.mat'),'pval'); 
            PlaceCells = find(pval > 0.95);
            n = length(PlaceCells); 
            
            for s=1:nCompSessions
                [r(:,:,s),~] = PlaceTimeCorr(mapMD,MD1,MDs(s),PlaceCells); 
            end
    end
        
    percentStable = squeeze(sum(r(:,2,:) < 0.05))./n;
end