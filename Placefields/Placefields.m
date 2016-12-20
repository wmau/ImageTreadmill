function Placefields(MD,varargin)
%Placefields(MD,varargin)
%
%   Main course function for calculating place fields, based off Dave
%   Sullivan's. Looks at epochs where mouse is running and produces a
%   spatial heatmap of firing for each neuron. Also runs a permutation test
%   to ask whether a place field is better than expected by chance. This is
%   done by shuffling transient times (deck of cards shuffle) while keeping
%   position fixed, misaligning position and spiking. Mutual information
%   (see Olypher et al., 2003 and spatInfo) is calculated using real and
%   shuffled data and compared to produce p-value.
%
%   INPUTS
%       MD: session to analyze. 
%
%       optional...
%       exclude_frames: vector of frames to exclude from aligned position
%       data. Can either be the indices of the aligned data OR the indices
%       of non-aligned TENASPIS data. Only tested this first case, but the
%       code for the second case should work.
%
%       cmperbin: centimeters per spatial bin. Default = 1.
%
%       minspeed: velocity threshold for saying that the mouse is running.
%       Default = 3 cm/s.
%
%       B: number of permutation iterations for determining the legitimacy
%       of a place field. Default = 10,000.
%
%       aligned: logical telling this function whether the Pos_data
%       variable you entered has already been aligned or not. Default =
%       true. 
%
%       Pos_data: output of PreprocessMousePosition_auto, aligned or not
%       aligned. Default = Pos_align.mat.
%
%       Tenaspis_data: output of Tenaspis. Default = FinalOutput.mat.

%% Parse inputs.
    cd(MD.Location);
    
    ip = inputParser;
    ip.addRequired('MD',@(x) isstruct(x)); 
    ip.addParameter('exclude_frames',[],@(x) isnumeric(x)); 
    ip.addParameter('cmperbin',1,@(x) isscalar(x)); 
    ip.addParameter('minspeed',3,@(x) isscalar(x)); 
    ip.addParameter('B',1000,@(x) isscalar(x));
    ip.addParameter('aligned',true,@(x) islogical(x));
    ip.addParameter('Pos_data','Pos_align.mat',@(x) ischar(x));
    ip.addParameter('Tenaspis_data','FinalOutput.mat',@(x) ischar(x)); 
    
    ip.parse(MD,varargin{:});
    
    %Compile.
    exclude_frames = ip.Results.exclude_frames;
    cmperbin = ip.Results.cmperbin;
    minspeed = ip.Results.minspeed;
    B = ip.Results.B; 
    aligned = ip.Results.aligned;
    Pos_data = ip.Results.Pos_data;
    Tenaspis_data = ip.Results.Tenaspis_data;
    
%% Set up.
    if aligned
        load('Pos_align.mat',...
            'PSAbool','x_adj_cm','y_adj_cm','speed','xmin','xmax','ymin','ymax'); 
        x = x_adj_cm; y = y_adj_cm; clear x_adj_cm y_adj_cm;
    else
        load('Pos.mat','xpos_interp','ypos_interp');
        load(Tenaspis_data,'PSAbool'); 
        [x,y,speed,PSAbool] = AlignImagingToTracking(MD.Pix2CM,PSAbool,0);
        xmin = min(x); ymin = min(y); 
        xmax = max(x); ymax = max(y);
        
        %Assuming your exclude_frames did not already apply to the aligned
        %data, correct them. This should work, but haven't actually tested
        %this.
        exclude_frames = exclude_frames - (FToffset-1);
        exclude_frames(exclude_frames < 0) = [];
        exclude_frames(exclude_frames > size(PSAbool,2)) = [];
    end
    
    %Basic variables. 
    [nNeurons,nFrames] = size(PSAbool); 
    velocity = convtrim(speed,ones(1,2*20))./(2*20);    %Smooth velocity (cm/s).
    good = true(1,nFrames);                             %Frames that are not excluded.
    good(exclude_frames) = false;
    isrunning = good;                                   %Running frames that were not excluded. 
    isrunning(velocity < minspeed) = false;
    
%% Get occupancy map. 
    lims = [xmin xmax;
            ymin ymax];
    [OccMap,RunOccMap,xEdges,yEdges,xBin,yBin] = ...
        MakeOccMap(x,y,lims,good,isrunning,cmperbin);

    %Don't need non-isrunning epochs anymore. 
    x = x(isrunning);
    y = y(isrunning);
    PSAbool = logical(PSAbool(:,isrunning));
    nGood = length(x); 
    
%% Construct place field and compute mutual information.
    %Preallocate.
    TCounts = cell(1,nNeurons);
    TMap_gauss = cell(1,nNeurons); 
    TMap_unsmoothed = cell(1,nNeurons); 
    pos = [x;y];
    parfor n=1:nNeurons    
        %Make place field.
        [TMap_unsmoothed{n},TCounts{n},TMap_gauss{n}] = ...
            MakePlacefield(PSAbool(n,:),pos,xEdges,yEdges,RunOccMap,...
            'cmperbin',cmperbin,'smooth',true);
    end
    
    %Compute mutual information.
    MI = spatInfo(TMap_unsmoothed,RunOccMap,PSAbool,true);
    
%% Get statistical significance of place field using mutual information.
    %Preallocate. 
    pval = nan(1,nNeurons);
    
    %Set up progress bar.
    resolution = 2;
    updateInc = round(nNeurons/(100/resolution));
    p = ProgressBar(100/resolution);
    parfor n=1:nNeurons
        
        %Predetermine transient frame shifts, disassociates transients from
        %location. 
        rTMap = cell(1,B);
        shifts = randi([0 nGood],B,1); 
        for i=1:B
            %Circular shift. 
            shuffled = circshift(PSAbool(n,:),[0 shifts(i)]);
            
            %Make place field from shifted transient vector. 
            rTMap{i} = MakePlacefield(shuffled,pos,xEdges,yEdges,...
                RunOccMap,'cmperbin',cmperbin,'smooth',false); 

        end

        %Calculate mutual information of randomized vectors. 
        rMI = spatInfo(rTMap,RunOccMap,repmat(PSAbool(n,:),[B,1]),false); 

        %Get p-value. 
        pval(n) = 1-(sum(MI(n)>rMI)/B); 
        
        if round(n/updateInc) == (n/updateInc)
            p.progress;
        end
    end
    p.stop; 
    
    save('Placefields.mat','OccMap','RunOccMap','TCounts','TMap_gauss',...
        'TMap_unsmoothed','minspeed','isrunning','cmperbin','exclude_frames',...
        'xEdges','yEdges','xBin','yBin','pval'); 
end