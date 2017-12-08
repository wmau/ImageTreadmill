function OnTreadmillMovie(md,clim,movietype,varargin)
%OnTreadmillMovie(animal,date,session,clim,movietype)
%       
%   Makes two synced AVIs. ImagingMovie_movietype is the calcium imaging
%   movie during treadmill run epochs. Blue outlines mark active neurons,
%   red outlines mark actiev time cells. TrackingMovie_movietype is the
%   concurrent Cineplex movie. Red dot indicates mouse position. 
%
%   INPUTS
%       animal, date, session: Example - 'GCamp6f_45_treadmill',
%       '11_20_2015', 2. 
%
%       clim: Color axis limits need to be specified beacause brightness
%       levels differ between mice. [0 2000] for G45 D1 movies, [0 4000]
%       for G48 D1 movies. 
%
%       movietype: string, either 'd1','smoothed','dff', or 'slpdf'. 
%

%% Preliminary stuff. 
    cd(md.Location);
    
    HalfWindow = 0; 
    neuraldata = 'FinalOutput.mat';
    load('MovieDims.mat','Xdim','Ydim');
    if ~isempty(varargin)
        if any(strcmp('halfwindow',varargin))   %HalfWindow (should be 10 for comparing D1 movies to T2 trace). 
            HalfWindow = varargin{find(strcmp('halfwindow',varargin))+1}; 
        end
        
        if any(strcmp('alt_input',varargin))    %For T2 outputs. 
            neuraldata = varargin{find(strcmp('alt_input',varargin))+1}; 
             
        end
        
        if any(strcmp('noi',varargin))
            highlight = varargin{find(strcmp('noi',varargin))+1};
        end
    end
    
    movietype = lower(movietype); 
    switch movietype
        case 'd1'
            h5file = fullfile(pwd,'D1Movie.h5');
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
    load(fullfile(pwd,neuraldata),'PSAbool','NeuronImage'); 
    nNeurons = length(NeuronImage); 
    %load(fullfile(pwd,'CC.mat'),'cc');
    load(fullfile(pwd,'TimeCells.mat'),'TodayTreadmillLog','TimeCells','movies'); 
    outlines = cellfun(@bwboundaries,NeuronImage,'unif',0); 
    colors = rand(nNeurons,3); 
    %Neuron centroids. 
    %centroids = getNeuronCentroids(animal,date,session); 
        
    %Align and get indices where mouse was on treadmill. 
    Pix2CM = 0.1256; sf = 0.6246;
    [x,y,~,~,FToffset,~,~,~] = AlignImagingToTracking(Pix2CM,PSAbool,HalfWindow); 
    x = x./Pix2CM*sf; 
    y = y./Pix2CM*sf; 
    treadmillInds = TodayTreadmillLog.inds;
    nRuns = size(treadmillInds,1); 
  
%% Initialize movie properties.     
    %Tracking movie. 
    trackingread = dir('*.avi'); 
    trackingread = trackingread.name; 
       
    %Imaging movie.
    imagingmovie = VideoWriter(['ImagingMovie_',movietype],'MPEG-4');
    imagingmovie.FrameRate = 20; 
    
    %Tracking movie.
    trackingread = VideoReader(trackingread); 
    trackingwrite = VideoWriter(['TrackingMovie_',movietype],'MPEG-4');
    trackingwrite.FrameRate = 20; 
    
    open(trackingwrite); 
    open(imagingmovie);
    
    tInc = 0;
    ifigure = figure('Position',[260 240 560 420]); 
    tfigure = figure('Position',[840 240 560 420]);
    for thisEpoch=[22 24 27 31 35]
        if TodayTreadmillLog.complete(thisEpoch)
            sFrame = treadmillInds(thisEpoch,1) + FToffset;
            eFrame = treadmillInds(thisEpoch,2) + FToffset; 

            for i=sFrame:eFrame
%% Imaging movie. 
                %Get frame. 
                frame = h5read(h5file,'/Object',[1 1 i+HalfWindow 1],[Xdim Ydim 1 1]);

                %Active neurons.
                active = find(PSAbool(:,i));

                %Display the imaging movie frame. 
                set(0,'currentfigure',ifigure); 
                imagesc(frame); caxis(clim); colormap gray;
                annotation(gcf,'textbox',[0.6, 0.8, 0.2, 0.07],'String',...
                {['t = ',num2str(round(tInc,1)),' seconds']},...
                'Color','white','EdgeColor','white','tag','iannotation');
                annotation(gcf,'textbox',[0.2, 0.8, 0.1, 0.07],'String',...
                {['Lap ',num2str(thisEpoch)]},'Color','white',...
                'EdgeColor','white','tag','iannotation');
            
                %If there are active neurons according to FT, plot its
                %outline.
                if ~isempty(active)
                    hold on;
                    for neuron=active'
                        if exist('highlight','var')
                            if ismember(neuron,highlight)                            
                                plot(outlines{neuron}{1}(:,2),outlines{neuron}{1}(:,1),...
                                    'Color',colors(neuron,:),'linewidth',2);
                            end
                        else
                            patchline(outlines{neuron}{1}(:,2),outlines{neuron}{1}(:,1),...
                                'edgecolor','k','edgealpha',0.2,'linewidth',2);
                        end
                    end
                    hold off;        
                end
                
                %Get rid of axis marks. 
                axis equal; axis off; 

                %Write video to the brain imaging movie. 
                F = getframe(gcf); 
                writeVideo(imagingmovie,F); 
                delete(findall(gcf,'Tag','iannotation')); 

%% Tracking movie
                %AVI file reads using time, not frames. 
                trackingread.currentTime = movies.t(treadmillInds(thisEpoch,1)) + tInc; 
                
                %Get frame. 
                frame = readFrame(trackingread); 
                
                %Display frame. 
                set(0,'currentfigure',tfigure); 
                imagesc(flipud(frame)); 
                hold on;
                t = findclosest(trackingread.currentTime,movies.aviFrame); 
                plot(x(t),y(t),'r.'); hold off; 
                annotation(gcf,'textbox',[0.6, 0.8, 0.2, 0.07],'String',...
                {['t = ',num2str(round(tInc,1)),' seconds']},'Color','red',...
                'EdgeColor','red','tag','tannotation');
                annotation(gcf,'textbox',[0.2, 0.8, 0.1, 0.07],'String',...
                {['Lap ',num2str(thisEpoch)]},'Color','red',...
                'EdgeColor','red','tag','tannotation');
                axis off; 
                
                %Write video to the tracking movie. 
                F = getframe(gcf);
                writeVideo(trackingwrite,F);  
                delete(findall(gcf,'Tag','tannotation')); 
                
                %At the end of a lap, write a few blank frames. 
                if i==eFrame
                    clf; 
                    F = getframe(gcf); 
                    
                    for j=1:30
                        writeVideo(imagingmovie,F);
                        writeVideo(trackingwrite,F); 
                    end
                end
                
                tInc = tInc + 0.05;     %Assumes 20 Hz sampling rate. Advance the time by 0.05s.
            end  
            
        end           
        
        %Reset increment. 
        tInc = 0; 
    end
    
    close(gcf); 
    close(imagingmovie); 
    close(trackingwrite); 
end