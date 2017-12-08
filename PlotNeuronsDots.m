function PlotNeuronsDots(md,varargin)

%
%
%

%%
    cd(md.Location); 
    load('MovieDims.mat','Xdim','Ydim'); 
    load('FinalOutput.mat','NumNeurons'); 
    
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x)); 
    p.addParameter('neurons',1:NumNeurons,@(x) isnumeric(x)); 
    p.addParameter('dotSizes',20,@(x) isnumeric(x));
    p.addParameter('dotColors','k'); 
    p.addParameter('filled',false,@(x) islogical(x));
    p.addParameter('transparency',1,@(x) isnumeric(x));
    p.addParameter('transformation',[]); 
    
    p.parse(md,varargin{:}); 
    neurons = p.Results.neurons;
    dotSizes = p.Results.dotSizes; 
    dotColors = p.Results.dotColors;
    filled = p.Results.filled; 
    transparency = p.Results.transparency; 
    transformation = p.Results.transformation; 
    
%% Get the centroids.
    %If we're not doing a transformation proceed normally. 
    if isempty(transformation)
        %If we already calculated the centroids, proceed. 
        try 
            load('NeuronMasks.mat','centroids')

        %Otherwise, find cell centroids and save so we don't have to do this
        %again. 
        catch
            centroids = getNeuronCentroids(md);

            %If this file already exists, append the centroid variable to it.
            %Otherwise, create the file. 
            save('NeuronMasks.mat','centroids'); 
        end

        %If we were able to load a NeuronMasks file but it doesn't have the
        %centroids variable in it, get the centroids then append that variable
        %to the mat file. 
        if ~exist('centroids','var')
            centroids = getNeuronCentroids(md); 
            save('NeuronMasks.mat','centroids','-append'); 
        end
        
    %Otherwise, we have to compute the centroids from scratch.
    else
        load('FinalOutput.mat','NeuronImage');
        
        %Transform the mask. 
        NeuronImage = cellfun(@(x) imwarp_quick(x,transformation),...
            NeuronImage,'UniformOutput',false); 
        
        %Basically rewriting the getNeuronCentroids function. 
        props = cellfun(@(x) regionprops(x,'Centroid'),NeuronImage(neurons));
        temp = extractfield(props,'Centroid'); 
        centroids = nan(NumNeurons,2);
        centroids(neurons,:) = [temp(1:2:end)', temp(2:2:end)'];
    end
        
%% Plot neurons as scatter. 
    if filled
        scatter(centroids(neurons,1),centroids(neurons,2),...
            dotSizes,dotColors,'filled','markerfacealpha',transparency);
    else 
        scatter(centroids(neurons,1),centroids(neurons,2),...
            dotSizes,dotColors,'markerfacealpha',transparency);
    end
    xlim([0 Ydim]); ylim([0 Xdim]);
    axis equal; axis off; 
end