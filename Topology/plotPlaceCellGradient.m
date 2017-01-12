function plotPlaceCellGradient(md)
%PlotTimeCellGradient(md)
%
%   Plots place cell ROIs in colors based on their rank in the sequence.
%   Colder colors are earlier, warmer colors are later. 
%
%   INPUT
%       md: Session entry. 

%% Plot outlines. 
    cd(md.Location);
    DATA = CompileMultiSessionData(md,{'placecells'});
    PlaceCells = DATA.placecells{1};
    
    load('Pos_align.mat','PSAbool');
    nNeurons = size(PSAbool,1);
    nPCs = length(PlaceCells); 
    c = colormap(jet(nPCs)); 
    [~,~,~,order] = LinearizedPFs_treadmill(md);
    PlaceCells = PlaceCells(order);
    
    PlotNeurons(md,1:nNeurons,[.7 .7 .7],1);
    for n=1:nPCs
        PlotNeurons(md,PlaceCells(n),c(n,:),2);
        hold on;
    end
end