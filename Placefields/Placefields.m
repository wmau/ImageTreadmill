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
%       of a place field. Default = 1,000.
%
%       aligned: logical telling this function whether the Pos_data
%       variable you entered has already been aligned or not. Default =
%       true. 
%
%       Pos_data: output of PreprocessMousePosition_auto, aligned or not
%       aligned. Default = Pos_align.mat.
%
%       Tenaspis_data: output of Tenaspis. Default = FinalOutput.mat.
%
%       exclude_frames: frames to exclude, either in raw imaging movie
%       (e.g. dropped frames), or in aligned data aligned data (e.g. frames
%       corresponding to times mouse is off the maze). If aligned = false,
%       then use raw imaging movie frames, if aligned = true, then use
%       frame numbers from data in Pos_align.mat
%
%       name_append: string to append to Placefields file ->
%       Placefieldsname_append.mat

%% Parse inputs.
    currdir = cd;

    [dirstr, MD] = ChangeDirectory(MD.Animal, MD.Date, MD.Session); % Change Directory and fill in partial MD if used
    
    ip = inputParser;
    ip.addRequired('MD',@(x) isstruct(x)); 
    ip.addParameter('exclude_frames',[],@(x) isnumeric(x)); 
    ip.addParameter('cmperbin',2.5,@(x) isscalar(x)); 
    ip.addParameter('minspeed',3,@(x) isscalar(x)); 
    ip.addParameter('B',1000,@(x) isscalar(x));
    ip.addParameter('aligned',true,@(x) islogical(x));
    ip.addParameter('Pos_data','Pos_align.mat',@(x) ischar(x));
    ip.addParameter('Tenaspis_data','FinalOutput.mat',@(x) ischar(x)); 
    ip.addParameter('name_append','',@ischar);
    ip.addParameter('rescale', 1, @isnumeric); %option to rescale data after converting to cm.
    ip.KeepUnmatched = true;
    ip.addParameter('saveSI',true,@(x) islogical(x)); 
    
    ip.parse(MD,varargin{:});
    
    %Compile.
    exclude_frames = ip.Results.exclude_frames;
    cmperbin = ip.Results.cmperbin;
    minspeed = ip.Results.minspeed;
    B = ip.Results.B; 
    aligned = ip.Results.aligned;
    Pos_data = ip.Results.Pos_data;
    Tenaspis_data = ip.Results.Tenaspis_data;
    name_append = ip.Results.name_append;
    saveSI = ip.Results.saveSI;
    rescale = ip.Results.rescale;
    
%% Set up.

    if aligned
        load(Pos_data,'PSAbool','x_adj_cm','y_adj_cm','speed','time_interp',...
            'xmin','xmax','ymin','ymax','FToffset'); 
        x = x_adj_cm; y = y_adj_cm; clear x_adj_cm y_adj_cm;
        offset = FToffset;
        
        % NRK code here to apply additional scaling factor if desired...
    else
        load(Pos_data,'xpos_interp','ypos_interp');
        load(Tenaspis_data,'PSAbool'); 
        [x,y,speed,PSAbool,offset,~,~,time_interp] = ...
            AlignImagingToTracking(MD.Pix2CM,PSAbool,0,'basedir',MD.Location);
        xmin = min(x); ymin = min(y); 
        xmax = max(x); ymax = max(y);
        
    end
    if ~isempty(rescale) && rescale ~= 1
        disp(['Rescaling position/speed data by a factor of ' num2str(rescale)])
        xmin = xmin*rescale;
        xmax = xmax*rescale;
        ymin = ymin*rescale;
        ymax = ymax*rescale;
        x = x*rescale;
        y = y*rescale;
        speed = speed*rescale;
    end
    
    % Adjust frames to exclude by FToffset
    if ~isempty(exclude_frames) && ~aligned
        exclude_frames = exclude_frames - (offset-2);
        exclude_frames(exclude_frames < 1) = [];
        exclude_frames(exclude_frames > size(PSAbool,2)) = [];
    end
    
   
    %Basic variables. 
    [nNeurons,nFrames] = size(PSAbool); 
    velocity = convtrim(speed,ones(1,2*20))./(2*20);    %Smooth velocity (cm/s).
    good = true(1,nFrames);                             %Frames that are not excluded.
    good(exclude_frames) = false;
    isrunning = good;                                   %Running frames that were not excluded. 
    isrunning(velocity < minspeed) = false;
    
%     keyboard
% %% De-bugging spot for aligning exclude_frames
%     
%     figure
% %     frames_plot = 1:length(x);
%     plot(time_interp,x,'k',time_interp(~good), x(~good),'r*')
    
%% Get occupancy map. 
    lims = [xmin xmax;
            ymin ymax];
    [OccMap,RunOccMap,xEdges,yEdges,xBin,yBin] = ...
        MakeOccMap(x,y,lims,good,isrunning,cmperbin);

    %Don't need non-isrunning epochs anymore. 
    x = x(isrunning);
    y = y(isrunning);
    PSAbool = logical(PSAbool(:,isrunning));
    xBin = xBin(isrunning);
    yBin = yBin(isrunning);
    nGood = length(x); 
    
%% Construct place field and compute mutual information.
    %Preallocate.
    TCounts = cell(1,nNeurons);
    TMap_gauss = cell(1,nNeurons); 
    TMap_unsmoothed = cell(1,nNeurons); 
    pos = [x;y];
    %%
    parfor n=1:nNeurons    
        %Make place field.
        [TMap_unsmoothed{n},TCounts{n},TMap_gauss{n}] = ...
            MakePlacefield(PSAbool(n,:),pos,xEdges,yEdges,RunOccMap,...
            'cmperbin',cmperbin,'smooth',true);
    end
    
    %Compute mutual information.
    MI = spatInfo(TMap_unsmoothed, RunOccMap, PSAbool, true, 'name_append', ...
        name_append);
    
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
        rMI = spatInfo(rTMap,RunOccMap,repmat(PSAbool(n,:),[B,1]),saveSI); 

        %Get p-value. 
        pval(n) = 1-(sum(MI(n)>rMI)/B); 
        
        if round(n/updateInc) == (n/updateInc)
            p.progress;
        end
    end
    p.stop; 
    
    save(fullfile(dirstr,['Placefields' name_append '.mat']),'OccMap','RunOccMap','TCounts','TMap_gauss',...
        'TMap_unsmoothed','minspeed','isrunning','cmperbin','exclude_frames',...
        'xEdges','yEdges','xBin','yBin','pval','x','y','PSAbool','MI'); 
    
    cd(currdir)
end