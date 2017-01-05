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
    load('TimeCells.mat','TimeCells','T'); 
    load('TemporalInfo.mat','sig');
    TimeCells = intersect(find(sig),TimeCells);
    
    load('Pos_align.mat','PSAbool');
    nNeurons = size(PSAbool,1);
    nTCs = length(TimeCells); 
    c = colormap(jet(nTCs)); 
    [~,order] = PastalkovaPlot(md,T,false);
    TimeCells = TimeCells(order);
    
    PlotNeurons(md,1:nNeurons,[.7 .7 .7],1);
    for n=1:nTCs
        PlotNeurons(md,TimeCells(n),c(n,:),2);
        hold on;
    end
end