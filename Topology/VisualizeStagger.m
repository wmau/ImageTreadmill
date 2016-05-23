function VisualizeStagger(md,graphData,neuron,varargin)
%
%
%

%% 
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x)); 
    p.addRequired('graphData',@(x) isstruct(x));
    p.addRequired('neuron',@(x) isnumeric(x) && isscalar(x)); 
    p.addParameter('edgelist',find(graphData.A(:,neuron))',@(x) isnumeric(x));
    p.addParameter('plotcells',0,@(x) islogical(x)); 
    p.parse(md,graphData,neuron,varargin{:});
    
    el = p.Results.edgelist; 
    A = p.Results.graphData.A;
    Ap = p.Results.graphData.Ap;
    null = p.Results.graphData.null;
    lagMat = p.Results.graphData.lagMat;
    plotcells = p.Results.plotcells;
    md = p.Results.md; 
    neuron = p.Results.neuron;
    
%%
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
    
    leadRaster = buildRaster(nLaps,inds,FT,neuron);
    
    %Line formats for raster. 
    lead = struct;
    lead.Color = 'r';
    lead.LineStyle = '-';
    
    lag = struct;
    lag.Color = 'g';
    lag.LineStyle = '-';

%% Plot neurons.
    nTicks = 6;
    if plotcells
        for i=el
            f = figure('Position',[185 70 1260 700]);
            %Neuron lagging. 
            subplot(3,5,1);
            imagesc([0:T],...
                [1:5:sum(TodayTreadmillLog.delaysetting==T)],...
                ratebylap(:,:,i)); 
            colormap gray; title(['Neuron ',num2str(i)]);
            xlabel('Time [s]'); ylabel('Laps');

            %Neuron leading. 
            subplot(3,5,2); 
            imagesc([0:T],...
                [1:5:sum(TodayTreadmillLog.delaysetting==T)],...
                ratebylap(:,:,neuron)); 
            colormap gray; title(['Neuron ',num2str(neuron)]);
            xlabel('Time [s]');  

            lagRaster = buildRaster(nLaps,inds,FT,i);

            %Raster. 
            subplot(3,5,6:7);
            plotSpikeRaster(lagRaster,'PlotType','vertline',...
                'LineFormat',lag,'TimePerBin',0.05,'SpikeDuration',0.05); 
            hold on;
            plotSpikeRaster(leadRaster,'PlotType','vertline',...
                'LineFormat',lead,'TimePerBin',0.05,'SpikeDuration',0.05); 
            ax = gca; 
            ax.XTick = linspace(ax.XLim(1),ax.XLim(2),nTicks);
            ax.XTickLabel = linspace(0,T,nTicks);
            ax.YTick = [1:5:nLaps];
            hold off; ylabel('Laps'); 

            %Histogram. 
            subplot(3,5,11:12); 
            histogram(null{i,neuron},'binwidth',1,'normalization','probability','facecolor','c'); 
            hold on;
            histogram(lagMat{i,neuron},'binwidth',1,'normalization','probability','facecolor','y'); 
            hold off;
            title(['P = ',num2str(Ap(i,neuron))]);
            xlabel('Time [s]'); ylabel('Proportion');

            %Topology
            subplot(3,5,[3:5,8:10,14:15]);
            PlotNeurons(md,[1:NumNeurons],'k',1);
            hold on;
            PlotNeurons(md,neuron,'r',2);
            PlotNeurons(md,i,'g',2);
            hold off
                  
            set(f,'PaperOrientation','landscape');
            set(f,'PaperUnits','normalized');
            set(f,'PaperPosition',[0 0 1 1]);
        end

%% Don't plot neurons. 
    else
        for i=el
            figure('Position',[185 70 490 700]);
            %Neuron lagging. 
            subplot(3,2,1);
            imagesc([0:T],...
                [1:5:sum(TodayTreadmillLog.delaysetting==T)],...
                ratebylap(:,:,i)); 
            colormap gray; title(['Neuron ',num2str(i)]);
            xlabel('Time [s]'); ylabel('Laps');

            %Neuron leading. 
            subplot(3,2,2); 
            imagesc([0:T],...
                [1:5:sum(TodayTreadmillLog.delaysetting==T)],...
                ratebylap(:,:,neuron)); 
            colormap gray; title(['Neuron ',num2str(neuron)]);
            xlabel('Time [s]');  

            lagRaster = buildRaster(nLaps,inds,FT,i);

            %Raster. 
            subplot(3,2,3:4);
            plotSpikeRaster(lagRaster,'PlotType','vertline',...
                'LineFormat',lag,'TimePerBin',0.05,'SpikeDuration',0.05); 
            hold on;
            plotSpikeRaster(leadRaster,'PlotType','vertline',...
                'LineFormat',lead,'TimePerBin',0.05,'SpikeDuration',0.05); 
            ax = gca; 
            ax.XTick = linspace(ax.XLim(1),ax.XLim(2),nTicks);
            ax.XTickLabel = linspace(0,T,nTicks);
            ax.YTick = [1:5:nLaps];
            hold off; ylabel('Laps'); 

            %Histogram. 
            subplot(3,2,5:6); 
            histogram(null{i,neuron},'binwidth',1,'normalization','probability','facecolor','c'); 
            hold on;
            histogram(lagMat{i,neuron},'binwidth',1,'normalization','probability','facecolor','y'); 
            hold off;
            title(['P = ',num2str(Ap(i,neuron))]);
            xlabel('Time [s]'); ylabel('Proportion');
        end
    end
end

function raster = buildRaster(nTrials,inds,FT,neuron) 
%
%
%

%%
    raster = zeros(nTrials,unique(diff(inds,[],2))+1);
    for t=1:nTrials
        lapRaster = [0 FT(neuron,inds(t,1):inds(t,2))]; 
        raster(t,:) = diff(lapRaster);
    end
    raster(raster == -1) = 0;
    raster = logical(raster);
end