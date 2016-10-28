function PlotNeurons(md,neurons,col,thickness)
%PlotNeurons(md,neurons,col,thickness)
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
    cd(md.Location); load('FinalOutput.mat','NeuronImage');
    
    w = whos('-file','FinalOutput.mat');
    if any(strcmp({w.name},'xOutline'))
        load('FinalOutput.mat','xOutline','yOutline'); 
    elseif exist('PlaceMaps.mat','file')
        load('PlaceMaps.mat','xOutline','yOutline'); 
    else
        NumNeurons = length(NeuronImage);

        xOutline = cell(NumNeurons,1);
        yOutline = cell(NumNeurons,1); 

        for i=1:NumNeurons
            b = bwboundaries(NeuronImage{i},'noholes');
            xOutline{i} = b{1}(:,1); 
            yOutline{i} = b{1}(:,2);
        end
        
        save('FinalOutput.mat','xOutline','yOutline','-append');
    end
    
    [Xdim,Ydim] = size(NeuronImage{1});
        
    hold on; 
    for i=neurons
        plot(yOutline{i},xOutline{i},'color',col,'linewidth',thickness); 
    end
    hold off; 
    xlim([0 Ydim]); ylim([0 Xdim]); 
    %line([0 100/1.16],[0 0],'linewidth',5,color,'k');
    axis equal; axis off;
     
end