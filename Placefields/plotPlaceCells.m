function plotPlaceCells(md,varargin)
%
%
%

%% Parse inputs.
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x));
    p.addParameter('neurons',[],@(x) isnumeric(x));
    
    p.parse(md,varargin{:});
    neurons = p.Results.neurons; 
    
%% Set up.
    cd(md.Location);
    
    if isempty(neurons)
        neurons = getPlaceCells(md,.01);
    end
    nPCs = length(neurons);

    load('Placefields.mat','TMap_gauss','isrunning');
    load('SpatialInfo.mat','MI','Ipos','okpix');
    load('PlacefieldStats.mat','PFarea','bestPF');
    
    PFarea = nansum(PFarea,2);
    try
        load('Pos_align.mat','PSAbool','x_adj_cm','y_adj_cm');
    catch
        load('FinalOutput.mat','PSAbool');
        [x_adj_cm,y_adj_cm,~,PSAbool] = AlignImagingToTracking(md.Pix2CM,PSAbool,0);
    end
    
    %For dotplot.
    PSAbool = logical(PSAbool);
    PSAbool(:,~isrunning) = false;
    
    %Rotate.
    [x,y] = rotate_arena(x_adj_cm,y_adj_cm,-90);
    clear x_adj_cm y_adj_cm;
    
    %Metrics.
    nNeurons = size(PSAbool,1);
    dims = size(TMap_gauss{1});
    
%% Flatten all the Ipos vectors. 
    IMap = cell(nNeurons,1);
    for n=1:nNeurons
        IMap{n} = unflattenIpos(Ipos{n},okpix,dims);
    end

%%  
    folder = 'C:\Users\William Mau\Documents\Projects\Time Cell Imaging Summer 2015 -\Paper\Figures\Supplementals\Example Cells\Place Cells';
    thisNeuron = 1;
    keepgoing = true;
    
    while keepgoing
        spks = PSAbool(neurons(thisNeuron),:);
        f = figure(50); 
        subplot(1,3,1);
            plot(x,y,'Color',[.7 .7 .7]); hold on;
            plot(x(spks),y(spks),'r.','MarkerSize',8);
            axis equal; axis off; 
            title(['Neuron #',num2str(neurons(thisNeuron))]);
            
        subplot(1,3,2);
            tmap = imagesc(TMap_gauss{neurons(thisNeuron)});
            set(tmap,'alphadata',~isnan(TMap_gauss{neurons(thisNeuron)}));
            axis equal; axis off; 
            colormap hot;         
            title([num2str(round(MI(neurons(thisNeuron)),2)),' bits']);
            
        subplot(1,3,3);
            imap = imagesc(IMap{neurons(thisNeuron)});
            set(imap,'alphadata',~isnan(TMap_gauss{neurons(thisNeuron)}));
            axis equal; axis off; 
            colormap hot;           
            title(['PF Area = ',num2str(PFarea(neurons(thisNeuron)))]);

        print(fullfile(folder,[md.Animal,' PC #',num2str(neurons(thisNeuron))]),'-dpdf');
            
        figure(50);
        [keepgoing,thisNeuron] = scroll(thisNeuron,nPCs,f);
        close all;
    end
end