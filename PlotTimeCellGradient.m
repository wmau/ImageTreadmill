function plotTimeCellGradient(md)
%PlotTimeCellGradient(md)
%
%   Plots time cell ROIs in colors based on their rank in the sequence.
%   Colder colors are earlier, warmer colors are later. 
%
%   INPUT
%       md: Session entry. 

%% Plot outlines. 
    cd(md.Location);
    load('TimeCells.mat','TimeCells'); 
    nTCs = length(TimeCells); 
    c = colormap(jet(nTCs)); 
    
    for n=1:nTCs
        PlotNeurons(md,TimeCells(n),c(n,:),1);
        hold on;
    end
end