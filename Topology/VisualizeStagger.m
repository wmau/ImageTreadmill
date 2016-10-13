function [triggerRaster,targetRaster,cellOffsetSpread,el] = VisualizeStagger(graphData,neuron,varargin)
%
%
%

%% Grab inputs. 
    p = inputParser;
    p.addRequired('graphData',@(x) isstruct(x));
    p.addRequired('neuron',@(x) isnumeric(x) && isscalar(x)); 
    p.addParameter('edgelist',find(graphData.A(:,neuron))',@(x) isnumeric(x));
    p.addParameter('plotcells',false,@(x) islogical(x)); 
    p.parse(graphData,neuron,varargin{:});
    
    el = p.Results.edgelist; 
    Atpval = p.Results.graphData.Atpval;
    Atrlpval = round(p.Results.graphData.Atrlpval,3);
    timeShuffleNulls = p.Results.graphData.tNullLats;
    trialShuffleNulls = p.Results.graphData.trlNullLats;
    latencies = p.Results.graphData.latencies;
    plotcells = p.Results.plotcells;
    neuron = p.Results.neuron;
    
    if isfield(graphData,'prune_p')
        graphData.p = graphData.prune_p;
    end
    
    if size(el,1) > size(el,2)
        el = el'; 
    end
    
%% Set up.
    %Get MD entry. 
    md = findMDfromGraphData(graphData); 
    
    %Change directory and load initial variables. 
    cd(md.Location); 
    load('TimeCells.mat','ratebylap','T','TodayTreadmillLog'); 
    load('Pos_align.mat','time_interp','FT');
    NumNeurons = size(FT,1);
    delays = TodayTreadmillLog.delaysetting; 
    complete = TodayTreadmillLog.complete;
    
    %Get treadmill run indices. 
    inds = getTreadmillEpochs(TodayTreadmillLog,time_interp); 
    inds = inds(find(complete),:);  %Only completed runs. 
    inds(:,2) = inds(:,1) + 20*T-1; %Consistent length.   
    
    %Trim ratebylap. 
    ratebylap = ratebylap(delays==T & complete,:,:);
    ratebylap = ratebylap(:,~isnan(ratebylap(1,:,1)),:); 
    nLaps = size(ratebylap,1); 
    
    %Build raster for second neuron. 
    targetRaster = buildRaster(inds,FT,neuron);
    
    %Line format for second neuron raster.
    lead.Color = 'r';
    lead.LineWidth = 1.5;
    
    %Line format for first neuron raster, all spikes, transparent. 
    lag.Color = [0 1 0 0.3];    %Transparent green.
    lag.LineWidth = 1.5;
    
    %Line format for first neuron raster, only spikes immediately
    %preceding those of neuron two. 
    immediatelag.Color = 'g';
    immediatelag.LineWidth = 1.5;
   
    %Preallocate.
    i=1;
    nInitiators = length(el);
    cellOffsetSpread = zeros(1,nInitiators);
    ratio = zeros(1,nInitiators);
    triggerRaster = cell(1,nInitiators);
    edges = 0:0.2:10;
%% Plot neurons.
    nTicks = 6;
    for e=el
        %% Raster - Neuron lagging. 
        f = figure('Position',[185 70 490 700]);
        subplot(3,2,1); 
        imagesc([0:T],...
            [1:nLaps],...
            ratebylap(:,:,e)); 
        colormap gray; title(['\color{green}Trigger \color{black}ROI #',num2str(e)]);
        ylabel('Laps');

        %% Raster - Neuron leading. 
        subplot(3,2,2);
        imagesc([0:T],...
            [1:nLaps],...
            ratebylap(:,:,neuron)); 
        colormap gray; title(['\color{red}Target \color{black}ROI #',num2str(neuron)]);
         
        %% Tick raster
        %Build the tick raster for neuron 1. 
        triggerRaster{i} = buildRaster(inds,FT,e);
        
        %Raster for responses immediately preceding neuron 2.
        [immediateRaster,d] = stripRaster(triggerRaster{i},targetRaster); 
       
        %Raster. 
        subplot(3,2,3:4);
        hold on;
        plotSpikeRaster(triggerRaster{i},'PlotType','vertline',...
            'LineFormat',lag,'VertSpikeHeight',1.5); 
        plotSpikeRaster(immediateRaster,'PlotType','vertline',...
            'LineFormat',immediatelag,'VertSpikeHeight',1.5); 
        plotSpikeRaster(targetRaster,'PlotType','vertline',...
            'LineFormat',lead,'VertSpikeHeight',1.5); 
        ax = gca; 
        ax.Color = 'k';
        ax.XTick = linspace(ax.XLim(1),ax.XLim(2),nTicks);
        ax.XTickLabel = linspace(0,T,nTicks);
        ax.YTick = [1:5:nLaps];
        set(gca,'ticklength',[0 0]);
        hold off; ylabel('Laps'); xlabel('Time [s]'); title('Paired Raster');

        %% Temporal distance histogram. 
        if plotcells, subplot(3,2,5); else subplot(3,2,5:6); end
        hold on;
                        
        %Time shuffle distribution.
        [timeNull,timeNullBins] = hist(timeShuffleNulls{e,neuron},edges); 
        timeNull = timeNull./length(timeShuffleNulls{e,neuron});
        stairs(timeNullBins,timeNull,'-.','linewidth',2,'color','c'); 
        
        %Trial shuffle distribution.
        [trialNull,trialNullBins] = hist(trialShuffleNulls{e,neuron},edges);
        trialNull = trialNull./length(trialShuffleNulls{e,neuron});
        stairs(trialNullBins,trialNull,'-.','linewidth',2,'color','b');
        
        %Treadmill latency.
        TMAlignedOnsets = TMLatencies(immediateRaster,targetRaster);
        [tmlat,tmlatBins] = hist(TMAlignedOnsets,edges);
        tmlat = tmlat./length(TMAlignedOnsets);
        stairs(tmlatBins,tmlat,':','linewidth',2,'color',[.6 .6 .6]);
        
        %MAD ratio.
        ratio(i) = round(mad(d,1)/ mad(TMAlignedOnsets,1),2);
        
        %Real latency distribution.
        [emp,empBins] = hist(latencies{e,neuron},edges);
        emp = emp./length(latencies{e,neuron});
        stairs(empBins,emp,'linewidth',2,'color','k'); 
        
        %Labels.
        title('Spike Time Latencies');
        xlabel('Latency from Target [s]'); ylabel('Proportion of Latencies');
        legend({['Time shuffle p = ',num2str(Atpval(e,neuron))],...
            ['Trial shuffle p = ',num2str(Atrlpval(e,neuron))],...
            ['Trdmll-trgt MAD ratio = ',num2str(ratio(i))],...
            'Trggr-trgt'},'fontsize',6);
        set(gca,'linewidth',1.5);

        %% Anatomical topology.
        if plotcells
            subplot(3,2,6);
            PlotNeurons(md,[1:NumNeurons],'k',1);
            hold on;
            PlotNeurons(md,neuron,'r',2);
            PlotNeurons(md,e,'g',2);
            hold off
                  
            %Sizing purpose for saving onto pdf. 
             set(f, 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPositionMode', 'auto');
        end
        
        %Advance counter.
        i = i+1; 
    end
    
end

%% stripRaster
