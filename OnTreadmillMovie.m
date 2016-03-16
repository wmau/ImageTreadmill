function OnTreadmillMovie(animal,date,session,clim,movietype,varargin)
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
%       clim: Color axis limits need to be specified beacause birghtness
%       levels differ between mice. [0 2000] for G45 D1 movies, [0 4000]
%       for G48 D1 movies. 
%
%       movietype: string, either 'D1' or 'smoothed', indicating whether to
%       use the first derivative movie or ICmovie_smoothed. 
%

%% Preliminary stuff. 
    close all;
    
    %h5 file. 
    ChangeDirectory(animal,date,session);
    
    HalfWindow = 0; 
    neuraldata = 'ProcOut.mat';
    if ~isempty(varargin)
        if any(strcmp('HalfWindow',varargin))   %HalfWindow (should be 10 for comparing D1 movies to T2 trace). 
            HalfWindow = find(strcmp('HalfWindow'),varargin))+1; 
        end
        
        if any(strcmp('alt_input',varargin))    %For T2 outputs. 
            neuraldata = find(strcmp('alt_input',varargin))+1; 
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
    end

%% Load data and process. 
    %Imaging data. 
    load(fullfile(pwd,neuraldata),'FT','NeuronImage'); 
    %load(fullfile(pwd,'CC.mat'),'cc');
    load(fullfile(pwd,'TimeCells.mat'),'TodayTreadmillLog','TimeCells','movies'); 
    outlines = cellfun(@bwboundaries,NeuronImage,'unif',0); 
    
    %Neuron centroids. 
    %centroids = getNeuronCentroids(animal,date,session); 
        
    %Align and get indices where mouse was on treadmill. 
    Pix2CM = 0.15; sf = 0.6246;
    [~,~,~,~,FToffset,~,~,~] = AlignImagingToTracking(Pix2CM,FT); 
    x = movies.x./Pix2CM*sf; y = movies.y./Pix2CM*sf; 
    treadmillInds = getTreadmillEpochs(TodayTreadmillLog,movies.t); 
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
    for thisEpoch=1:nRuns
        if TodayTreadmillLog.complete(thisEpoch)
            sFrame = treadmillInds(thisEpoch,1) + FToffset;
            eFrame = treadmillInds(thisEpoch,2) + FToffset; 

            for i=sFrame:eFrame
%% Imaging movie. 
                %Get frame. 
                try 
                    frame = h5read(h5file,'/Object',[1 1 i+HalfWindow 1],[Xdim Ydim 1 1]);
                catch
                    disp([movietype,' movie not found! Try another type.']); 
                end

                %Active neurons.
                active = find(FT(:,i));

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
                        if ismember(neuron,TimeCells), c = '-r';
                        else c = '-b'; end       
                        plot(outlines{neuron}{1}(:,2),outlines{neuron}{1}(:,1),c,'linewidth',2);
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
                t = findclosest(trackingread.currentTime,movies.t); 
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