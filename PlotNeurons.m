function PlotNeurons(md,neurons,col,thickness,varargin)
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

%% 
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x));
    p.addRequired('neurons',@(x) isnumeric(x)); 
    p.addRequired('col',@(x) isnumeric(x) | ischar(x)); 
    p.addRequired('thickness',@(x) isnumeric(x)); 
    p.addParameter('transformation',[]); 
    p.addParameter('lineStyle','-',@(x) ischar(x)); 
    p.addParameter('transparent',false,@(x) islogical(x)); 
    
    p.parse(md,neurons,col,thickness,varargin{:}); 
    transformation = p.Results.transformation; 
    lineStyle = p.Results.lineStyle; 
    transparent = p.Results.transparent; 
    
%% Plot neurons. 
    cd(md.Location); 
    load('MovieDims.mat','Xdim','Ydim');
    load('FinalOutput.mat','NumNeurons');
    
    if size(neurons,1) > size(neurons,2)
        neurons = neurons';
    end
        
    try
        load('NeuronMasks.mat','xOutline','yOutline');
    catch
        load('FinalOutput.mat','NeuronImage');
        
        if ~isempty(transformation)
            NeuronImage = cellfun(@(x) imwarp_quick(x,transformation),...
                NeuronImage,'UniformOutput',false); 
        end

        %Preallocate.
        xOutline = cell(NumNeurons,1);
        yOutline = cell(NumNeurons,1); 

        %Make the outlines of the cell ROIs. 
        for i=1:NumNeurons
            b = bwboundaries(NeuronImage{i},'noholes');
            try
                xOutline{i} = b{1}(:,1); 
                yOutline{i} = b{1}(:,2);

            %Handles the rare case where a transformation takes a cell out
            %of the field of view. 
            catch
                xOutline{i} = nan;
                yOutline{i} = nan;
            end
        end

        %Only save the non-transformed layout. We don't want to keep
        %translating cells. 
        if isempty(transformation) 
            if exist('NeuronMasks.mat','file')~=2
                save('NeuronMasks.mat','xOutline','yOutline');
            else
                save('NeuronMasks.mat','xOutline','yOutline','-append'); 
            end
        end
    end
    
    if ~isempty(transformation)
        load('FinalOutput.mat','NeuronImage');

        NeuronImage = cellfun(@(x) imwarp_quick(x,transformation),...
            NeuronImage,'UniformOutput',false); 
        
        %Preallocate.
        xOutline = cell(NumNeurons,1);
        yOutline = cell(NumNeurons,1); 

        %Make the outlines of the cell ROIs. 
        for i=1:NumNeurons
            b = bwboundaries(NeuronImage{i},'noholes');
            try
                xOutline{i} = b{1}(:,1); 
                yOutline{i} = b{1}(:,2);

            %Handles the rare case where a transformation takes a cell out
            %of the field of view. 
            catch
                xOutline{i} = nan;
                yOutline{i} = nan;
            end
        end
    end
    
    %If only one row of RGB values is detected, replicate it. 
    if size(col,1) == 1
        col = repmat(col,NumNeurons,1); 
    end
    
    hold on; 
    
    for i=neurons
        cellOutline = plot(yOutline{i},xOutline{i},'color',col(i,:),'linewidth',thickness,...
            'linestyle',lineStyle); 
        
        if transparent
            cellOutline.Color(4) = 0.3; 
        end
    end
    hold off; 
    xlim([0 Ydim]); ylim([0 Xdim]); 
    %line([0 100*1.10],[0 0],'linewidth',5,'color','k');
    axis equal; axis off;
     
end