function TCStills(md,varargin)
%
%
%

%% Set up.
    %Get some things initialized. 
    cd(md.Location); 
    load('TimeCells.mat','T','TodayTreadmillLog','TimeCells'); 
    load('FinalOutput.mat','xOutline','yOutline','NeuronImage');
    load('Pos_align.mat','FToffset');
    [inds,nRuns] = TrimTrdmllInds(TodayTreadmillLog,T);
    totalNeurons = length(xOutline);
    load('MovieDims.mat','Xdim','Ydim');
    
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x));
    p.addParameter('neurons',TimeCells,@(x) isnumeric(x)); 
    p.addParameter('edges',0:1:10,@(x) isnumeric(x));
    p.addParameter('movie','D1Movie.h5',@(x) ischar(x)); 
    p.addParameter('wndw',20,@(x) isscalar(x)); 
    p.addParameter('laps',1:nRuns,@(x) isnumeric(x));
    p.addParameter('plotit',true,@(x) islogical(x) | ismember(x,[0 1]));
    
    p.parse(md,varargin{:}); 
    neurons = p.Results.neurons;
    nNeurons = length(neurons);
    edges = p.Results.edges;
    movie = p.Results.movie; 
    edgeSize = length(edges)-1;
    wndw = p.Results.wndw;
    plotit = p.Results.plotit;
    
    if strcmp(movie,'D1movie.h5'), HalfWindow = 10; else HalfWindow = 0; end;
    
    centroids = getNeuronCentroids(md);
%% 
    %Set conditions still extraction time frames. 
    resolution = max(edges)/(length(edges)-1);
    chunkSize = resolution * 20; 

    %Set up loop.
    stills = cell(1,totalNeurons);
    borderStarts = zeros(totalNeurons,2); 
    borderLims = zeros(totalNeurons,2); 
    pRes = 2; 
    updateInc = round(nNeurons/(100/pRes));  
    p = ProgressBar(100/pRes);
    for n = 1:nNeurons
        nn = neurons(n);        %Neuron #.
        
        %Border defining h5 movie window. There is some confusing stuff
        %going on here. In h5read we want the row and column indices which
        %line on the Y and X axes, respectively. That's why we flip the X,Y
        %coordinates for centroids here.
        borderStarts(nn,:) = round(fliplr(centroids(nn,:)) - wndw);
        borderLims(nn,:) = wndw*2*ones(1,2);

        %Fix edge cases. 
        exceed = borderStarts(nn,:) < 1;
        
        %Xdim is the number of ROWS in the frame, which is in the Y
        %dimension, corresponding to the first value of borderStarts, which
        %specifies ROWS. 
        if (borderStarts(nn,1) + borderLims(nn,1)) > Xdim
            borderLims(nn,1) = floor(Xdim - borderStarts(nn,1)); 
        elseif exceed(1)       
            borderLims(nn,1) = borderLims(nn,1) + borderStarts(nn,1);
        end 
        
        %Ydim, similarly, is the number of COLUMNS, which is in the X
        %dimension, corresponding to the second value of borderStarts,
        %which specifies COLUMNS.
        if (borderStarts(nn,2) + borderLims(nn,2)) > Ydim
            borderLims(nn,2) = floor(Ydim - borderStarts(nn,2));
        elseif exceed(2) 
            borderLims(nn,2) = borderLims(nn,2) + borderStarts(nn,2);
        end
        
        %Fix edge cases.
        borderStarts(nn,exceed) = 1; 
            
        %Set up loop for each time bin.
        stills{nn} = nan([borderLims(nn,:),edgeSize]);
        for t = 1:edgeSize
            
            %Preallocate an empty array for the frame average. 
            temp = zeros([borderLims(nn,:),nRuns]);
            for r = 1:nRuns
                %Time index.
                i = inds(r,1) + edges(t)*20 + HalfWindow + FToffset;
                
                %Read h5 file and take the mean of the chunk across time.
                temp(:,:,r) = mean(h5read(movie,'/Object',[borderStarts(nn,:) i 1],...
                     [borderLims(nn,:) chunkSize 1]),3);
            end

            %Take the mean across laps.
            stills{nn}(:,:,t) = nanmean(temp,3);
            
        end
        
        %Update progress bar. 
        if round(n/updateInc) == (n/updateInc)
            p.progress;
        end
    end
    p.stop;   

%%
    keepgoing = true;
    n = 1;
    if plotit
        while keepgoing
            thisNeuron = neurons(n);
            
            %Define the x and y borders. Remember x here is ROWS (so
            %actually the y dimension on imagesc). Conversely, y is
            %COLUMNS, so the x dimension on imagesc.
            x = borderStarts(thisNeuron,1):(borderStarts(thisNeuron,1)+(borderLims(thisNeuron,1)-1));
            y = borderStarts(thisNeuron,2):(borderStarts(thisNeuron,2)+(borderLims(thisNeuron,2)-1));
            
            %Make mask then find the color axis limits by looking for the
            %min and max within that mask.
            mask = logical(NeuronImage{thisNeuron}(x,y));       
            cLims = getMaxMin(stills{thisNeuron},mask,1:edgeSize);
            
            f = figure(thisNeuron);       
            for t=1:edgeSize
                %Plot.
                h = subplot_auto(edgeSize,t);
                hold on;
                imagesc(y,x,stills{thisNeuron}(:,:,t));                 %Frame.    
                set(gca,'YDir','reverse');                              %Flip imagesc default.
                plot(yOutline{thisNeuron},xOutline{thisNeuron},'r',...
                    'linewidth',2);                                     %ROI outline.
                caxis([cLims]); axis tight; axis off;
                text(y(2),x(end-5),[num2str(edges(t)),'-'...
                    num2str(edges(t+1)),'s'],'color','c','fontsize',14);
            end
            colormap gray;
            
            
            [keepgoing,n] = scroll(n,length(neurons),f);
            close all;
        end
    end
end

function cLims = getMaxMin(stills,mask,inds)
%cLims = getMaxMin(stills,mask,inds)
%
%  Gets the maximum and minimum of a three dimensional matrix masked by a
%  binary mask.
%
%   INPUTS
%       stills: XxYxT matrix of frames.
%
%       mask: neuron mask.
%
%       inds: indices of time.
%

%% 
    nInds = length(inds);
    maskedMaxMins = zeros(nInds,2); 
    for i = 1:nInds
        still = stills(:,:,inds(i)); 
        maskedMaxMins(i,1) = max(max(still(mask)));
        maskedMaxMins(i,2) = min(min(still(mask)));
    end
    
    cLims = [min(maskedMaxMins(:,1)) max(maskedMaxMins(:,2))];
end