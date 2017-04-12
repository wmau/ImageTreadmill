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
    
    %Color for plotting. 
    purple = [.58 .44 .86];
%% Match cells.
    map = msMatchCells(getMapMD(mds),mds,neurons,false); 

%% Get all rasters and tuning curves. 
    DATA = CompileMultiSessionData(mds,{'si','pfcorr'});
    SI = DATA.si; 
    PFCORR = DATA.pfcorr; 
    
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
            %Get neuron ID. 
            neuron = map(thisNeuron,thisSession);
            
            %Raster. 
            rasterAX(thisSession) = subplot(nSessions,2,thisSession*2-1);
                imagesc(rasters{thisSession}(:,:,thisNeuron)); 
                title(['Neuron #',num2str(neuron)]);
                colormap hot; caxis([0 1.2]);
                set(gca,'xtick',[],'ytick',[1,size(rasters{thisSession},1)]);
                ylabel('Laps'); 

            %Tuning curve. 
            curveAX(thisSession) = subplot(nSessions,2,thisSession*2);
                plot(bins,smoothCurve{thisSession}(thisNeuron,:),...
                    'color',purple,'linewidth',2);              
                ylabel('Rate');     
                xlim([0 max(bins)]);
                set(gca,'linewidth',4,'xtick',[0 max(bins)],'xticklabel',...
                    {'0','163'});
                
            %Label date or whether neuron was detected. Also label spatial
            %information.
            if isnan(neuron) || neuron==0
                title('Not detected');       
                ylim([0 eps]);
            else
                %Label day. 
                title(['Day ',num2str(thisSession)]);               
            end
            if thisSession==nSessions, xlabel('Linearized Distance (cm)'); end 
                
        end
        
        %Normalize axes to max among all days. 
        curveYLims = [0, max([curveAX.YLim])];
        set(curveAX,'YLim',curveYLims);
        
        %Label statistics. 
        for thisSession=1:nSessions
            %Get neuron ID. 
            neuron = map(thisNeuron,thisSession);
            neuronIsGood = ~isnan(neuron) && neuron~=0;
            
            %Go to subplot then label spatial information. 
            subplot(nSessions,2,thisSession*2);
            
            if neuronIsGood
                text(0.2*max(bins),0.9*curveYLims(2),['SI = ',...
                    num2str(round(SI{thisSession}(neuron),3)),' bits']); 
            end

            %Label correlation coefficient. 
            if thisSession~=nSessions && neuronIsGood
                %Get critical p-value. 
                PCs = getPlaceCells(mds(thisSession),0.01);
                PCs = EliminateUncertainMatches([mds(thisSession),mds(thisSession+1)],PCs);

                crit = 0.01/length(PCs);
                if PFCORR{thisSession}(neuron,2) < crit
                    c = purple; 
                else, c = 'r'; 
                end

                %Label correlation coefficient. 
                text(0.6*max(bins),-.25*curveYLims(2),...
                    ['R = ',num2str(round(PFCORR{thisSession}(neuron,1),3))],...
                    'color',c); 
            end
        end

        %Scroll through neurons. 
        [keepgoing,thisNeuron] = scroll(thisNeuron,nNeurons,f);
        close all;
        
    end
end