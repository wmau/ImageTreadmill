function [triggerRaster,targetRaster,cellOffsetSpread,el] = VisualizeStagger_alt(md,graphData,neuron,direction,varargin)
%[triggerRaster,targetRaster,cellOffsetSpread,el] = VisualizeStagger_alt(md,graphData,neuron,direction,varargin)
%
%

%% Grab inputs. 
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x)); 
    p.addRequired('graphData',@(x) isstruct(x));
    p.addRequired('neuron',@(x) isnumeric(x) && isscalar(x)); 
    p.addRequired('direction',@(x) ischar(x)); 
    
    if strcmp(direction,'left'), lr = 1; else lr = 2; end;
    
    p.addParameter('edgelist',find(graphData.A{lr}(:,neuron))',@(x) isnumeric(x));
    p.addParameter('plotcells',false,@(x) islogical(x)); 
    p.parse(md,graphData,neuron,direction,varargin{:});
    
    el = p.Results.edgelist; 
    Ap = p.Results.graphData.Ap{lr};
    nulld = p.Results.graphData.nulld{lr};
    CC = p.Results.graphData.CC{lr};
    plotcells = p.Results.plotcells;
    md = p.Results.md; 
    neuron = p.Results.neuron;
    
%% Set up.
    %Change directory and load initial variables. 
    cd(md.Location); 
    load('TimeCells.mat','ratebylap','T','TodayTreadmillLog'); 
    load('Pos_align.mat','aviFrame','FT');
    load('Alternation.mat'); 
    goodLaps = Alt.summary(:,2) == lr & Alt.summary(:,3);
    NumNeurons = size(FT,1);
    delays = TodayTreadmillLog.delaysetting; 
    complete = TodayTreadmillLog.complete;
    
    %Get treadmill run indices. 
    inds = getTreadmillEpochs(TodayTreadmillLog,aviFrame); 
    inds = inds(find(delays==T & complete),:);  %Only completed runs. 
    inds(:,2) = inds(:,1) + 20*T-1; %Consistent length.   
    
    %Sanity check - are all the trial numbers for complete treadmill runs
    %unique? 
    trialsOnTM = Alt.trial(inds(:,1)); 
    if length(unique(trialsOnTM)) ~= length(trialsOnTM)
        disp('Warning! Multiple laps sorted by postrials_treadmill correspond to run indices!'); 
        
        trialCounts = hist(trialsOnTM,unique(trialsOnTM)); 
        badTrials = find(trialCounts>1); 
        inds(badTrials,:) = [];
    end
    
    %Trim ratebylap. 
    ratebylap = ratebylap(delays==T & complete,:,:);
    ratebylap(badTrials,:,:) = [];
    ratebylap = ratebylap(goodLaps,:,:);
    ratebylap = ratebylap(:,~isnan(ratebylap(1,:,1)),:); 
    nLaps = size(ratebylap,1); 
    
    %Build raster for second neuron. 
    targetRaster = buildRaster(inds,FT,neuron);
    targetRaster = targetRaster(goodLaps,:);
    
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
    triggerRaster = cell(1,nInitiators);
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
        colormap gray; title(['\color{green}Trigger \color{black}ROI #',num2str(e)]);
        ylabel('Laps');

        %% Raster - Neuron leading. 
        if plotcells, subplot(3,5,2); else subplot(3,2,2); end
        imagesc([0:T],...
            [1:nLaps],...
            ratebylap(:,:,neuron)); 
        colormap gray; title(['\color{red}Target \color{black}ROI #',num2str(neuron)]);
         

        %% Tick raster
        %Build the tick raster for neuron 1. 
        triggerRaster{i} = buildRaster(inds,FT,e);
        triggerRaster{i} = triggerRaster{i}(goodLaps,:); 
        
        %Raster for responses immediately preceding neuron 2.
        [immediateRaster,d] = stripRaster(triggerRaster{i},targetRaster); 
       
        %Raster. 
        if plotcells, subplot(3,5,6:7); else subplot(3,2,3:4); end
        plotSpikeRaster(triggerRaster{i},'PlotType','vertline',...
            'LineFormat',lag,'TimePerBin',0.05,'SpikeDuration',0.05); 
        hold on;
        plotSpikeRaster(immediateRaster,'PlotType','vertline',...
            'LineFormat',immediatelag,'TimePerBin',0.05,'SpikeDuration',0.05); 
        plotSpikeRaster(targetRaster,'PlotType','vertline',...
            'LineFormat',lead,'TimePerBin',0.05,'SpikeDuration',0.05); 
        ax = gca; 
        ax.Color = 'k';
        ax.XTick = linspace(ax.XLim(1),ax.XLim(2),nTicks);
        ax.XTickLabel = linspace(0,T,nTicks);
        ax.YTick = [1:5:nLaps];
        set(gca,'ticklength',[0 0]);
        hold off; ylabel('Laps'); xlabel('Time [s]'); 

        %% Temporal distance histogram. 
        if plotcells, subplot(3,5,11); else subplot(3,2,5); end
        histogram(-nulld{e,neuron},[0:0.25:10],'normalization','probability',...
            'facecolor','c'); 
        hold on;
        histogram(-CC{e,neuron},[0:0.25:10],'normalization','probability',...
            'facecolor','y'); 
        hold off;
        title({'Spike Time Latencies',...
            ['P = ',num2str(Ap(e,neuron))]});
        xlabel('Latency from Target [s]'); ylabel('Proportion of Spike Pairs');
        legend({'Shuffled','Trigger'});
        set(gca,'linewidth',1.5);

        %% Activity relative to cell vs relative to treadmill.
        %Only look at laps where both neurons were active. Immediate raster
        %only has trues on laps where leadRaster was active. 
        TMAlignedOnsets = TMLatencies(immediateRaster,targetRaster);

        %Spread of responses relative to treadmill start. 
        treadmillOffsetSpread = mad(TMAlignedOnsets,1);
        
        if plotcells, subplot(3,5,12); else subplot(3,2,6); end
        histogram(TMAlignedOnsets,[0:0.25:10],'normalization','probability',...
            'facecolor','k');
        hold on;      
        
        %Get spread. 
        cellOffsetSpread(i) = mad(d,1);
        
        %Ratio between cell-to-cell vs cell-to-treadmill.
        ratio(i) =  cellOffsetSpread(i) / treadmillOffsetSpread;
        
        %Histogram.
        histogram(-d,[0:0.25:10],'normalization','probability',...
            'facecolor','y');
        title({'Trigger-Target vs. Treadmill-Target',...
            ['TT Score = ',num2str(ratio(i))]});
        legend({'Treadmill','Trigger'});
        xlabel('Latency from Target [s]')
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
