function msLinearizedPF_raster(mds,varargin)
%msLinearizedPF_raster(mds,varargin)
%
%   Plot linearized rasters and tuning curves of spatial responses in place
%   cells.
%

%% Parse inputs.
    p = inputParser;
    p.addRequired('mds',@(x) isstruct(x)); 
    p.addParameter('minspeed',3,@(x) isscalar(x));
    p.addParameter('nBins',80,@(x) isscalar(x)); 
    p.addParameter('neurons',[],@(x) isnumeric(x)); 
    
    p.parse(mds,varargin{:}); 
    
    minspeed = p.Results.minspeed; 
    nBins = p.Results.nBins; 
    neurons = p.Results.neurons; 
    
    if isempty(neurons)
        neurons = getPlaceCells(mds(1),.01);
    end
    nNeurons = length(neurons);
    nSessions = length(mds);
    
%% Match cells.
    map = msMatchCells(getMapMD(mds),mds,neurons,false); 

%% Get all rasters and tuning curves. 
    [rasters,smoothCurve] = deal(cell(nSessions,1));
    for s=1:nSessions
        [rasters{s},smoothCurve{s}] = LinearizedPF_raster(mds(s),'neurons',map(:,s),...
            'plotit',false,'minspeed',minspeed); 
    end

%% Plot.
    thisNeuron = 1;
    keepgoing = true;
    bins = [1:0.001:nBins]';
    while keepgoing
        %Make the figure. 
        f = figure(40);
        f.Position = [-1300 -40 370 180*nSessions];
        
        %Plot each session.
        for thisSession=1:nSessions
            %Raster. 
            rasterAX(thisSession) = subplot(nSessions,2,thisSession*2-1);
                imagesc(rasters{thisSession}(:,:,thisNeuron)); 
                colormap hot; caxis([0 1.2]);
                set(gca,'xtick',[],'ytick',[1,size(rasters{thisSession},1)]);
                ylabel('Laps'); 

            %Tuning curve. 
            curveAX(thisSession) = subplot(nSessions,2,thisSession*2);
                plot(bins,smoothCurve{thisSession}(thisNeuron,:),'color',[.58 .44 .86],'linewidth',5);              
                ylabel('Rate');
                yLims = get(gca,'ylim');
                ylim([0 yLims(2)]);
                xlim([0 max(bins)]);
                set(gca,'linewidth',4);
                
            %Label. 
            if isnan(map(thisNeuron,thisSession)) || map(thisNeuron,thisSession)==0
                title('Not detected'); 
            else
                title(['Day ',num2str(thisSession)]);
            end
            if thisSession==nSessions, xlabel('Linearized Distance'); end 
        end
        
        %Normalize axes to max among all days. 
        curveXLims = [min([curveAX.XLim]), max([curveAX.XLim])];
        curveYLims = [min([curveAX.YLim]), max([curveAX.YLim])];
        set(curveAX,'XLim',curveXLims,'YLim',curveYLims);
            
        %Scroll through neurons. 
        [keepgoing,thisNeuron] = scroll(thisNeuron,nNeurons,f);
    end
end