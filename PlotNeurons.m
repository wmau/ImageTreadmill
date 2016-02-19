function PlotNeurons(NeuronImage,Xdim,Ydim,c,thickness)
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
    nNeurons = length(NeuronImage);
    
    hold on; 
    for i=1:nNeurons
        b = bwboundaries(NeuronImage{i});
        
        x = b{1}(:,1);
        y = b{1}(:,2);
        
        plot(y,x,'b','color',c,'linewidth',thickness); 
    end
    hold off; 
    xlim([0 Xdim]); ylim([0 Ydim]); axis off; 
     
end