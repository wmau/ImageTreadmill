function msVisualizeStagger(mapMD,md,neuron,A,varargin)
%msVisualizeStagger(mapMD,md,neuron,A,varargin)
%
%   Plots rasters of the neuron specified as the input plus all its
%   triggers. Also plots the latency histogram of the target from its
%   triggers and the treadmill. Then it looks for the same neuron across
%   the days specified in md and does the same thing. 
%
%   INPUTS
%       mapMD: session entry where the neuron map lives. 
%
%       md: session entries. The first entry must be the one that
%       corresponds to the other inputs A and neuron. Panels in the final
%       plot will appear in the same order as specified here. 
%
%       neuron: target neuron you wish to examine. 
%
%       A: adjacency matrix, from graphData.
%
%       varargin: 
%           -edgelist: list of triggers (or not; specifying ROIs will force
%           this function to plot them).
%

%% Setup. 
    %Parse inputs. 
    p = inputParser;
    p.addRequired('mapMD',@(x) isstruct(x));
    p.addRequired('md',@(x) isstruct(x));
    p.addRequired('neuron',@(x) isnumeric(x) && isscalar(x));
    p.addRequired('A',@(x) islogical(x) || isnumeric(x));
    p.addParameter('edgelist',find(A(:,neuron))',@(x) isnumeric(x));
    
    p.parse(mapMD,md,neuron,A,varargin{:});
    el = p.Results.edgelist;
    nTriggers = length(el);
    nSessions = length(md);
       
    %Get data from all the sessions then find the corresponding cells. 
    DATA = CompileMultiSessionData(md,{'ratebylap','ft','ttl','t','A'});
    targets = msMatchCells(mapMD,md,neuron);
    triggers = msMatchCells(mapMD,md,el);
    
    %Line format for second neuron raster.
    lead.Color = 'r';
    lead.LineWidth = 2;
    
    %Line format for first neuron raster, all spikes, transparent. 
    lag.Color = [0 1 0 0.3];    %Transparent green.
    lag.LineWidth = 2;
    
    %Line format for first neuron raster, only spikes immediately
    %preceding those of neuron two. 
    immediatelag.Color = 'g';
    immediatelag.LineWidth = 2;
    
    nTicks = 6;
     
    for e=1:nTriggers
        f = figure('Position',[80 240 380*nSessions 570]);
        
        goodSessions = find(targets > 0 & triggers(e,:) > 0);
        
        for s=goodSessions
            nLaps = sum(DATA.ttl{s}.complete);
            
            %% Raster - Trigger
            subplot(3,2*nSessions,s*2-1);
            imagesc([0:DATA.t{s}],...
                [1:nLaps],...
                DATA.ratebylap{s}(:,:,triggers(e,s))); 
            colormap gray; title(['\color{green}Trigger \color{black}ROI #',num2str(triggers(e,s))]);
            if s==1, ylabel('Laps'); end

            %% Raster - target
            subplot(3,2*nSessions,s*2);
            imagesc([0:DATA.t{s}],...
                [1:nLaps],...
                DATA.ratebylap{s}(:,:,targets(s)));
            colormap gray; title(['\color{red}Target \color{black}ROI #',num2str(targets(s))]);

            %% Tick raster
            inds = DATA.ttl{s}.inds(logical(DATA.ttl{s}.complete));
            inds(:,2) = inds(:,1) + 20*DATA.t{s}-1;
            triggerRaster = buildRaster(inds,DATA.ft{s},triggers(e,s));
            targetRaster = buildRaster(inds,DATA.ft{s},targets(s));

            [immediateRaster,d] = stripRaster(triggerRaster,targetRaster);

            subplot(3,2*nSessions,s*2+nSessions*2-1:s*2+nSessions*2);
            plotSpikeRaster(triggerRaster,'PlotType','vertline',...
                'LineFormat',lag,'VertSpikeHeight',1.5); 
            hold on;
            plotSpikeRaster(immediateRaster,'PlotType','vertline',...
                'LineFormat',immediatelag,'VertSpikeHeight',1.5); 
            plotSpikeRaster(targetRaster,'PlotType','vertline',...
                'LineFormat',lead,'VertSpikeHeight',1.5);   
                ax = gca; 
                ax.Color = 'k';
                ax.XTick = linspace(ax.XLim(1),ax.XLim(2),nTicks);
                ax.XTickLabel = linspace(0,DATA.t{s},nTicks);
                ax.YTick = [1:5:nLaps];
                set(gca,'ticklength',[0 0]);
                hold off; if s==1, ylabel('Laps'); end

            %% Target to trigger/treadmill latency
            %Get treadmill-target latencies. 
            TMAlignedOnsets = TMLatencies(immediateRaster,targetRaster); 

            %Spread of responses relative to treadmill start. 
            treadmillOffsetSpread = mad(TMAlignedOnsets,1);

            subplot(3,2*nSessions,s*2+nSessions*4-1:s*2+nSessions*4);
            histogram(TMAlignedOnsets,[0:0.25:10],'normalization','probability',...
                'facecolor','k');
            hold on;      

            %Get spread. 
            cellOffsetSpread = mad(d,1);

            %Ratio between cell-to-cell vs cell-to-treadmill.
            ratio =  cellOffsetSpread / treadmillOffsetSpread;

            %Histogram.
            histogram(-d,[0:0.25:10],'normalization','probability',...
                'facecolor','y');
            title(['TT = ',num2str(ratio)]);
            if DATA.A{s}(triggers(s),targets(s))
                title(['\color{green}TT = ',num2str(ratio)]);
            end
            
            legend({'Treadmill','Trigger'});
            xlabel('Latency from Target [s]'); 
            if s==1, ylabel('Proportion of Spikes'); end
            set(gca,'linewidth',1.5);
            hold off;
        end
            
        set(f,'PaperOrientation','landscape');
%         set(f,'PaperUnits','normalized');
%         set(f,'PaperPosition',[0 0 1 1]);
    end
        
end