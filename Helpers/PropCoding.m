function proportion = PropCoding(md,cellType)
%
%
%

%%
    cd(md.Location); 
    load('FinalOutput.mat','NumNeurons'); 
    
    codingCells = AcquireTimePlaceCells(md,cellType); 
    proportion = length(codingCells)/NumNeurons; 
end