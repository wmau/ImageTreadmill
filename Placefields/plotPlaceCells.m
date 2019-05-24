function plotPlaceCells(md,varargin)
%
%
%

%% Parse inputs.

    
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x));
    p.addParameter('neurons',[],@(x) isnumeric(x));
    p.addParameter('crit', 0.01, @(a) a > 0 && a <= 1);
    p.addParameter('nHits', 10, @(a) round(a,0) == a);
    p.addParameter('ratioHits', 0.2, @(a) a > 0 && a <= 1);
    p.addParameter('name_append','',@ischar);
    
    
    p.parse(md,varargin{:});
    neurons = p.Results.neurons; 
    crit = p.Results.crit;
    name_append = p.Results.name_append;
    
%% Set up.
    [dirstr, md] = ChangeDirectory(md.Animal, md.Date, md.Session); % Change Directory and fill in partial MD if used
    
    if isempty(neurons)
        neurons = getPlaceCells(md, crit, varargin{:});
    end
    nPCs = length(neurons);
    
    if nPCs == 0
        error('No neurons meet your placefield criteria!')
    end


    try
        load(fullfile(dirstr,['Placefields' name_append '.mat']),'TMap_gauss','isrunning');
    catch
        Placefields(md,'name_append', name_append);
    end
    load(fullfile(dirstr,['SpatialInfo' name_append '.mat']),'MI','Ipos','okpix');
    load(fullfile(dirstr,['PlacefieldStats' name_append '.mat']),'PFarea')
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
        try
        IMap{n} = unflattenIpos(Ipos{n},okpix,dims);
        catch
            keyboard
        end
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

        %print(fullfile(folder,[md.Animal,' PC #',num2str(neurons(thisNeuron))]),'-dpdf');
                    
        [keepgoing,thisNeuron] = scroll(thisNeuron,nPCs,f);
        close all;
    end
end