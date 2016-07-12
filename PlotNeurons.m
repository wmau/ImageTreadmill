function PlotNeurons(md,neurons,col,thickness)
%PlotNeurons(NeuronImage)
%
%   Given cell arrays of binary images, plot the outlines of neurons. 
%
%   INPUT
%       NeuronImage: Output from Tenaspis. 
%
%   OUTPUT
%       Plot of neuron outlines. 
%

%% Plot neurons. 
    cd(md.Location); 
    load('FinalOutput.mat','NeuronImage');
    [Xdim,Ydim] = size(NeuronImage{1});
    NumNeurons = length(NeuronImage);
    
    w = whos('file','FinalOutput.mat');
    if any(strcmp({w.name},'xOutline'))
        load('FinalOutput.mat','xOutline','yOutline'); 
    elseif exist('PlaceMaps.mat','file')
        load('PlaceMaps.mat','xOutline','yOutline'); 
    else
        xOutline = cell(NumNeurons,1);
        yOutline = cell(NumNeurons,1); 

        for i=1:NumNeurons
            b = bwboundaries(NeuronImage{i});
            xOutline{i} = b{1}(:,1); 
            yOutline{i} = b{1}(:,2);
        end
        
        save('FinalOutput.mat','xOutline','yOutline','-append');
    end
        
    hold on; 
    for i=neurons
        plot(yOutline{i},xOutline{i},'color',col,'linewidth',thickness); 
    end
    hold off; 
    %xlim([0 Xdim]); ylim([0 Ydim]); 
    axis off; 
     
end