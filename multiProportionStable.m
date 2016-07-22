function [percentStable,nStable] = multiProportionStable(mapMD,MD1,MDs,neurontype)
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
    load(fullfile(MD1.Location,'FinalOutput.mat'),'FT'); 
    NumNeurons = size(FT,1); 
    neurontype = lower(neurontype); 
    r = zeros(NumNeurons,2,nCompSessions); 
    switch neurontype
        case 'time'
            load(fullfile(MD1.Location,'TimeCells.mat'),'TimeCells'); 
            n = length(TimeCells);
            
            for s=1:nCompSessions
                [~,r(:,:,s)] = PlaceTimeCorr(mapMD,MD1,MDs(s),TimeCells); 
                
                %[~,~,~,MAP,MAPcols] = nMatchedNeurons(mapMD,MD1,MDs(s));
                %n(s) = sum(ismember(MAP(:,MAPcols(2)),TimeCells));
                n(s) = sum(~isnan(r(:,2,s)));
            end
        case 'place'
            load(fullfile(MD1.Location,'PlaceMaps.mat'),'pval'); 
            PlaceCells = find(pval > 0.95);
            n = length(PlaceCells); 
            
            for s=1:nCompSessions
                [r(:,:,s),~] = PlaceTimeCorr(mapMD,MD1,MDs(s),PlaceCells);
                
                %[~,~,~,MAP,MAPcols] = nMatchedNeurons(mapMD,MD1,MDs(s));
                %n(s) = sum(ismember(MAP(:,MAPcols(2)),PlaceCells));
                n(s) = sum(~isnan(r(:,2,s)));
            end
    end
        
    n=n';
    nStable = squeeze(sum(r(:,2,:) < 0.05));
    percentStable = nStable./n;
end