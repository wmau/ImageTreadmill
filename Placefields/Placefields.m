function Placefields(MD,varargin)
%
%
%

%% Parse inputs.
    cd(MD.Location);
    
    p = inputParser;
    p.addRequired('MD',@(x) isstruct(x)); 
    p.addParameter('exclude_frames',[],@(x) isnumeric(x)); 
    p.addParameter('cmperbin',1,@(x) isscalar(x)); 
    p.addParameter('minspeed',3,@(x) isscalar(x)); 
    p.addParameter('HalfWindow',0,@(x) isscalar(x)); 
    p.addParameter('B',500,@(x) isscalar(x));
    p.addParameter('Tenaspis_output','FinalOutput.mat',@(x) ischar(x)); 
    
    p.parse(MD,varargin{:});
    
    %Compile.
    exclude_frames = p.Results.exclude_frames;
    cmperbin = p.Results.cmperbin;
    minspeed = p.Results.minspeed;
    HalfWindow = p.Results.HalfWindow; 
    B = p.Results.B; 
    Tenaspis_output = p.Results.Tenaspis_output;
    
%% 
    load('Pos_align.mat','FT','x_adj_cm','y_adj_cm','speed','xmin','xmax','ymin','ymax'); 
    x = x_adj_cm; y = y_adj_cm; clear x_adj_cm y_adj_cm;
    load(Tenaspis_output,'NeuronImage','NeuronPixels'); 
    
    [nNeurons,nFrames] = size(FT); 
    velocity = convtrim(speed,ones(1,2*20))./(2*20);
    good = true(1,nFrames); 
    good(exclude_frames) = false;
    isrunning = good; 
    isrunning(velocity < minspeed) = false;
    
%% Get occupancy map. 
    lims = [xmin xmax;
            ymin ymax];
    [OccMap,RunOccMap,xEdges,yEdges] = MakeOccMap(x,y,lims,good,isrunning,cmperbin);

%% 
    TCounts = cell(1,nNeurons);
    TMap_gauss = cell(1,nNeurons); 
    TMap_unsmoothed = cell(1,nNeurons); 
    for n=1:nNeurons
        spkpos = [  FT(n,isrunning);...
                    x(isrunning);...
                    y(isrunning)];
        if sum(spkpos(1,:)) > 4
            [TCounts{n},TMap_gauss{n},TMap_unsmoothed{n}] = ...
                MakePlacefield(spkpos,xEdges,yEdges,RunOccMap,'cmperbin',cmperbin);
        end
    end
    
    save('Placefields.mat','OccMap','RunOccMap','TCounts','TMap_gauss',...
        'TMap_unsmoothed','minspeed','isrunning','cmperbin','exclude_frames'); 
end