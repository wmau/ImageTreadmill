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
    cd(md.Location); 
    load('MovieDims.mat','Xdim','Ydim');
    
    load('FinalOutput.mat','xOutline','yOutline');
    if ~exist('xOutline','var')
        load('FinalOutput.mat','NeuronImage');
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
        
    hold on; 
    for i=neurons
        plot(yOutline{i},xOutline{i},'color',col,'linewidth',thickness); 
    end
    hold off; 
    xlim([0 Ydim]); ylim([0 Xdim]); 
    %line([0 100*1.10],[0 0],'linewidth',5,color,'k');
    axis equal; axis off;
     
end