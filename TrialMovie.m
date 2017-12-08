function TrialMovie(md,clim,movietype,varargin)
%TrialMovie(animal,date,session,clim,movietype)
%       
%   Makes a movie of single trials. Imaging movie on the top left, tracking
%   movie on the bottom right. Also displays time from treadmill start and
%   lap #. Will highlight requested neurons in distinguishable colors if
%   provided. Otherwise, will highlight all active neurons in gray. Default
%   save name is TrialMovie_*movietype*.mp4.
%
%   INPUTS
%       md: session entry
%
%       clim: Color axis limits need to be specified beacause brightness
%       levels differ between mice. [0 2000] for G45 D1 movies, [0 4000]
%       for G48 D1 movies. 
%
%       movietype: string, either 'd1','smoothed','dff','slpdf', or
%       'bpdff'.
%
%       NAME,VALUE:
%           neurons: vector of neurons you want to plot. 
%
%           trials: vector of trials. 
%

%% Preliminary stuff. 
    cd(md.Location);
    
    p = inputParser; 
    p.addRequired('md',@(x) isstruct(x));
    p.addRequired('clim',@(x) isnumeric(x)); 
    p.addRequired('movietype',@(x) ischar(x)); 
    p.addParameter('neurons',[],@(x) isnumeric(x));
    p.addParameter('trials',[],@(x) isnumeric(x));
    
    p.parse(md,clim,movietype,varargin{:});
    
    neurons = p.Results.neurons; 
    if size(neurons,1) > size(neurons,2), neurons = neurons'; end
    
    HalfWindow = 0;
    movietype = lower(movietype); 
    load('MovieDims','Xdim','Ydim');
    switch movietype
        case 'd1'
            h5file = fullfile(pwd,'D1Movie.h5');
            HalfWindow = 10;
        case 'smoothed'
            cd('ICmovie_smoothed-Objects'); 
            h5file = fullfile(pwd,'Obj_1 - ICmovie_smoothed.h5'); 
            cd ..
        case 'dff'
            h5file = fullfile(pwd,'DFF.h5'); 
        case 'slpdf'
            h5file = fullfile(pwd,'SLPDF.h5'); 
        case 'bpdff'
            h5file = fullfile(pwd,'BPDFF.h5');
    end

%% Load data and process. 
    %Imaging data. 
    load(fullfile(pwd,'Pos_align.mat'),'FToffset');
    load(fullfile(pwd,'FinalOutput.mat'),'NeuronImage','PSAbool');
    load(fullfile(pwd,'Alternation.mat')); 
    parsed = Alt; clear Alt; 
    
    nNeurons = length(NeuronImage); 
    %load(fullfile(pwd,'CC.mat'),'cc');
    load(fullfile(pwd,'TimeCells.mat'),'TodayTreadmillLog','TimeCells','movies'); 
    outlines = cellfun(@bwboundaries,NeuronImage,'unif',0); 
    colors = distinguishable_colors(nNeurons); 
    %Neuron centroids. 
    %centroids = getNeuronCentroids(animal,date,session); 
        
    %Align and get indices where mouse was on treadmill. 
    Pix2CM = 0.1256; sf = 0.6246;
    [x,y] = AlignImagingToTracking(Pix2CM,PSAbool,HalfWindow); 
    load(fullfile(pwd,'Pos_align.mat'),'PSAbool');
    x = x./Pix2CM*sf; 
    y = y./Pix2CM*sf; 
    
    %Get indices for each trial. 
    nTrials = max(parsed.summary(:,1)); 
    trialInds = nan(nTrials,2); 
    for r=1:nTrials
        trialInds(r,1) = find(parsed.trial==r,1);
        trialInds(r,2) = find(parsed.trial==r,1,'last'); 
    end
    treadmillInds = TodayTreadmillLog.inds; 
  
%% Initialize movie properties.     
    %Tracking movie. 
    trackingread = dir('*.avi'); 
    trackingread = trackingread.name; 
    trackingread = VideoReader(trackingread); 
       
    %Movie.
    movieFile = VideoWriter(['TrialMovie_',movietype],'MPEG-4');
    movieFile.FrameRate = 20; 
    open(movieFile); 
    
    tInc = 0;
    figure('Position',[390 20 1100 960]); 
    for thisEpoch=[22 24 27 31 35]
        if TodayTreadmillLog.complete(thisEpoch)
            sFrame = trialInds(thisEpoch,1);        %Start frame.
            eFrame = trialInds(thisEpoch,2);        %End frame. 

            for i=sFrame:eFrame
%% Imaging movie. 
                %Get frame. 
                frame = h5read(h5file,'/Object',[1 1 i+HalfWindow+FToffset 1],[Xdim Ydim 1 1]);

                %Active neurons.
                active = find(PSAbool(:,i));

                %Display the imaging movie frame. 
                subplot(3,3,[1,2,4,5]); 
                imagesc(frame); caxis(clim); colormap gray;
                
                hold on;
                for neuron=neurons
                    patchline(outlines{neuron}{1}(:,2),outlines{neuron}{1}(:,1),...
                        'edgecolor',colors(neuron,:),'edgealpha',.4,'linewidth',2);
                end
            
                %If there are active neurons according to FT, plot its
                %outline.
                if ~isempty(active)
                    hold on;
                    for neuron=active'
                        if ~isempty(neurons)
                            if ismember(neuron,neurons)                            
                                plot(outlines{neuron}{1}(:,2),outlines{neuron}{1}(:,1),...
                                    'Color',colors(neuron,:),'linewidth',6); 
                            end
                        %Otherwise, plot outlines of all neurons in light
                        %gray. 
                        else
                            patchline(outlines{neuron}{1}(:,2),outlines{neuron}{1}(:,1),...
                                'edgecolor','k','edgealpha',0.2,'linewidth',2);
                        end
                    end
                    hold off;        
                end
                
                %Get rid of axis marks. 
                axis equal; axis off;                

%% Tracking movie
                %AVI file reads using time, not frames. 
                currentTime = movies.t(trialInds(thisEpoch,1)) + tInc;  
                trackingread.currentTime = currentTime;
                
                timeReltoTrdmll = currentTime - (movies.t(treadmillInds(thisEpoch,1))); 
                %Get frame. 
                frame = readFrame(trackingread); 
                
                %Display frame. 
                subplot(3,3,9);
                imagesc(flipud(frame)); set(gca,'ydir','reverse');
                hold on;
                t = findclosest(trackingread.currentTime,movies.aviFrame); 
                plot(x(t),y(t),'r.'); hold off; 
                set(gca,'ydir','normal');
                annotation(gcf,'textbox',[0.7, 0.8, 0.3, 0.02],'String',...
                {['Lap ',num2str(thisEpoch)],['t = ',num2str(round(timeReltoTrdmll,1)),...
                ' seconds']},'Color','k','EdgeColor','none','tag','tannotation',...
                'fontsize',20);
                axis off; 
        
                yPos = 0.4; 
                inc = 0;
                for neuron=neurons
                    annotation(gcf,'textbox',[0.1 yPos-inc 0.3 0.02],'String',...
                        ['Neuron ',num2str(neuron)],'Color',colors(neuron,:),...
                        'EdgeColor','none','fontsize',20,'tag','tannotation');
                    inc = inc+0.03;
                end

                %Legend for figure in paper. Comment the above loop and
                %uncomment these lines to reproduce. 
%                 annotation(gcf,'textbox',[0.1 0.37 0.3 0.02],'String',...
%                     'Fig. 1d left','Color',colors(906,:),...
%                     'EdgeColor','none','fontsize',20,'tag','tannotation');
% 
%                 annotation(gcf,'textbox',[0.1 0.34 0.3 0.02],'String',...
%                     'Fig. 1d right','Color',colors(431,:),...
%                     'EdgeColor','none','fontsize',20,'tag','tannotation');
%                 
%                 annotation(gcf,'textbox',[0.1 0.31 0.3 0.02],'String',...
%                     'Supplementary Fig. 1c #1','Color',colors(240,:),...
%                     'EdgeColor','none','fontsize',20,'tag','tannotation');
%                 
%                 annotation(gcf,'textbox',[0.1 0.28 0.3 0.02],'String',...
%                     'Fig. 1f left','Color',colors(152,:),...
%                     'EdgeColor','none','fontsize',20,'tag','tannotation');
%                 
%                 annotation(gcf,'textbox',[0.1 0.25 0.3 0.02],'String',...
%                     'Fig. 1f right','Color',colors(190,:),...
%                     'EdgeColor','none','fontsize',20,'tag','tannotation');
%                 
%                 annotation(gcf,'textbox',[0.1 0.22 0.3 0.02],'String',...
%                     'Supplementary Fig. 1d #1','Color',colors(84,:),...
%                     'EdgeColor','none','fontsize',20,'tag','tannotation');
                
                 
                %Write video to the tracking movie. 
                F = getframe(gcf);
                writeVideo(movieFile,F);  
                delete(findall(gcf,'Tag','tannotation')); 
                
                %At the end of a lap, write a few blank frames. 
                if i==eFrame
                    clf; 
                    F = getframe(gcf); 
                    
                    for j=1:30
                        writeVideo(movieFile,F);
                    end
                end
                
                tInc = tInc + 0.05;     %Assumes 20 Hz sampling rate. Advance the time by 0.05s.
            end  
            
        end           
        
        %Reset increment. 
        tInc = 0; 
    end
    
    close(gcf); 
    close(movieFile); 
end