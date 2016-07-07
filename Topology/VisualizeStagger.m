function [lagRaster,leadRaster,cellOffsetSpread,el] = VisualizeStagger(md,graphData,neuron,varargin)
%
%
%

%% Grab inputs. 
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x)); 
    p.addRequired('graphData',@(x) isstruct(x));
    p.addRequired('neuron',@(x) isnumeric(x) && isscalar(x)); 
    p.addParameter('edgelist',find(graphData.A(:,neuron))',@(x) isnumeric(x));
    p.addParameter('plotcells',false,@(x) islogical(x)); 
    p.parse(md,graphData,neuron,varargin{:});
    
    el = p.Results.edgelist; 
    Ap = p.Results.graphData.Ap;
    null = p.Results.graphData.null;
    lagMat = p.Results.graphData.lagMat;
    closest = p.Results.graphData.closest;
    plotcells = p.Results.plotcells;
    md = p.Results.md; 
    neuron = p.Results.neuron;
    
%% Set up.
    %Change directory and load initial variables. 
    cd(md.Location); 
    load('TimeCells.mat','ratebylap','T','TodayTreadmillLog'); 
    load('Pos_align.mat','aviFrame','FT');
    load('ProcOut.mat','NumNeurons');
    delays = TodayTreadmillLog.delaysetting; 
    complete = TodayTreadmillLog.complete;
    
    %Get treadmill run indices. 
    inds = getTreadmillEpochs(TodayTreadmillLog,aviFrame); 
    inds = inds(find(complete),:);  %Only completed runs. 
    inds(:,2) = inds(:,1) + 20*T-1; %Consistent length.   
    
    %Trim ratebylap. 
    ratebylap = ratebylap(delays==T & complete,:,:);
    ratebylap = ratebylap(:,~isnan(ratebylap(1,:,1)),:); 
    nLaps = size(ratebylap,1); 
    
    %Build raster for second neuron. 
    leadRaster = buildRaster(inds,FT,neuron);
    
    %Line format for second neuron raster.
    lead.Color = 'r';
    lead.LineWidth = 1;
    
    %Line format for first neuron raster, all spikes, transparent. 
    lag.Color = [0 1 0 0.3];    %Transparent green.
    lag.LineWidth = 2;
    
    %Line format for first neuron raster, only spikes immediately
    %preceding those of neuron two. 
    immediatelag.Color = 'g';
    immediatelag.LineWidth = 2;
   
    %Preallocate.
    i=1;
    nInitiators = length(el);
    cellOffsetSpread = zeros(1,nInitiators);
    ratio = zeros(1,nInitiators);
    lagRaster = cell(1,nInitiators);
%% Plot neurons.
    nTicks = 6;
    for e=el
        %% Raster - Neuron lagging. 
        if plotcells, windowWidth = 1260; else windowWidth = 490; end
        f = figure('Position',[185 70 windowWidth 700]);
        if plotcells, subplot(3,5,1); else subplot(3,2,1); end
        imagesc([0:T],...
            [1:nLaps],...
            ratebylap(:,:,e)); 
        colormap gray; title(['Neuron ',num2str(e)]);
        ylabel('Laps');

        %% Raster - Neuron leading. 
        if plotcells, subplot(3,5,2); else subplot(3,2,2); end
        imagesc([0:T],...
            [1:nLaps],...
            ratebylap(:,:,neuron)); 
        colormap gray; title(['Neuron ',num2str(neuron)]);
         

        %% Tick raster
        %Build the tick raster for neuron 1. 
        lagRaster{i} = buildRaster(inds,FT,e);
        
        %Raster for responses immediately preceding neuron 2.
        [immediateRaster,d] = stripRaster(lagRaster{i},leadRaster); 
       
        %Raster. 
        if plotcells, subplot(3,5,6:7); else subplot(3,2,3:4); end
        plotSpikeRaster(lagRaster{i},'PlotType','vertline',...
            'LineFormat',lag,'TimePerBin',0.05,'SpikeDuration',0.05); 
        hold on;
        plotSpikeRaster(immediateRaster,'PlotType','vertline',...
            'LineFormat',immediatelag,'TimePerBin',0.05,'SpikeDuration',0.05); 
        plotSpikeRaster(leadRaster,'PlotType','vertline',...
            'LineFormat',lead,'TimePerBin',0.05,'SpikeDuration',0.05); 
        ax = gca; 
        ax.Color = 'k';
        ax.XTick = linspace(ax.XLim(1),ax.XLim(2),nTicks);
        ax.XTickLabel = linspace(0,T,nTicks);
        ax.YTick = [1:5:nLaps];
        hold off; ylabel('Laps'); xlabel('Time [s]'); 

        %% Temporal distance histogram. 
        if plotcells, subplot(3,5,11); else subplot(3,2,5); end
        histogram(null{e,neuron},[-10:0.5:10],'normalization','probability',...
            'facecolor','c'); 
        hold on;
        histogram(lagMat{e,neuron},[-10:0.5:10],'normalization','probability',...
            'facecolor','y'); 
        hold off;
        title(['P = ',num2str(Ap(e,neuron))]);
        xlabel('Time [s]'); ylabel('Proportion');
        set(gca,'linewidth',1.5);

        %% Activity relative to cell vs relative to treadmill.
        %Only look at laps where both neurons were active. Immediate raster
        %only has trues on laps where leadRaster was active. 
        bothActiveLaps = find(any(immediateRaster,2));  
        TMAlignedOnsets = [];
        for l=bothActiveLaps'
            %Get the onset times of each neuron. 
            TMAlignedOnsets = [TMAlignedOnsets find(leadRaster(l,:))];    
        end
        %[~,TMalignedOnsets] = find(leadRaster);

        %Divide by frame rate. 
        TMAlignedOnsets = TMAlignedOnsets./20;

        %Spread of responses relative to treadmill start. 
        treadmillOffsetSpread = mad(TMAlignedOnsets,1);
        
        if plotcells, subplot(3,5,12); else subplot(3,2,6); end
        histogram(TMAlignedOnsets,[0:0.5:10],'normalization','probability',...
            'facecolor','k');
        hold on;      
        
        %Get spread. 
        cellOffsetSpread(i) = mad(d,1);
        
        %Ratio between cell-to-cell vs cell-to-treadmill.
        ratio(i) = cellOffsetSpread(i) / treadmillOffsetSpread;
        
        %Histogram.
        histogram(-d,[0:0.5:10],'normalization','probability',...
            'facecolor','y');
        title({['Spread Ratio = ',num2str(ratio(i))],...
            ['N = ', num2str(length(d)), ' comparisons']});
        xlabel('Temporal Distance [s]')
        set(gca,'linewidth',1.5);
        
        %% Anatomical topology.
        if plotcells
            subplot(3,5,[3:5,8:10,14:15]);
            PlotNeurons(md,[1:NumNeurons],'k',1);
            hold on;
            PlotNeurons(md,neuron,'r',2);
            PlotNeurons(md,e,'g',2);
            hold off
                  
            %Sizing purpose for saving onto pdf. 
            set(f,'PaperOrientation','landscape');
            set(f,'PaperUnits','normalized');
            set(f,'PaperPosition',[0 0 1 1]);
        end
        
        %Advance counter.
        i = i+1; 
    end
    
end

%% stripRaster
