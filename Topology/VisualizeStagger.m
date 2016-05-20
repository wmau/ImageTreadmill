function VisualizeStagger(md,graphData,neuron)
%
%
%

%% 
    %Change directory and load initial variables. 
    cd(md.Location); 
    load('TimeCells.mat','ratebylap','T','TodayTreadmillLog'); 
    load('Pos_align.mat','aviFrame','FT');
    load('ProcOut.mat','NumNeurons');
    delays = TodayTreadmillLog.delaysetting; 
    complete = TodayTreadmillLog.complete;
    A = graphData.A; 
    Ap = graphData.Ap;
    null = graphData.null;
    lagMat = graphData.lagMat; 
    
    %Get treadmill run indices. 
    inds = getTreadmillEpochs(TodayTreadmillLog,aviFrame); 
    inds = inds(find(complete),:);  %Only completed runs. 
    inds(:,2) = inds(:,1) + 20*T-1; %Consistent length.   
    
    %Trim ratebylap. 
    ratebylap = ratebylap(delays==T & complete,:,:);
    ratebylap = ratebylap(:,~isnan(ratebylap(1,:,1)),:); 
    nLaps = size(ratebylap,1); 
    
    leadRaster = buildRaster(nLaps,inds,FT,neuron);
    
    %Get edge list.
    el = find(A(:,neuron)); 
    
    %Line formats for raster. 
    lead = struct;
    lead.Color = 'r';
    lead.LineStyle = '-';
    
    lag = struct;
    lag.Color = 'g';
    lag.LineStyle = '-';
    
    for i=el'
        figure('Position',[185 70 1260 700]);
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
        ax.XTickLabel = [0:1:T];
        ax.YTick = [1:5:nLaps];
        hold off; ylabel('Laps'); 
        
        %Histogram. 
        subplot(3,5,11:12); 
        histogram(null{i,neuron},'binwidth',1,'normalization','probability'); 
        hold on;
        histogram(lagMat{i,neuron},'binwidth',1,'normalization','probability'); 
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