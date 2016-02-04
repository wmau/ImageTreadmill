function OnTreadmillMovie(animal,date,session,clim,movietype)
%
%
%

%% Preliminary stuff. 
    close all;
    
    %h5 file. 
    ChangeDirectory(animal,date,session);
    
    movietype = lower(movietype); 
    switch movietype
        case 'd1'
            h5file = fullfile(pwd,'D1Movie.h5');
        case 'smoothed'
            cd('ICmovie_smoothed-Objects'); 
            h5file = fullfile(pwd,'Obj_1 - ICmovie_smoothed.h5'); 
            cd ..
    end

%% Load data and process. 
    %Imaging data. 
    load(fullfile(pwd,'ProcOut.mat'),'FT','NeuronImage','Xdim','Ydim'); 
    load(fullfile(pwd,'CC.mat'),'cc');
    load(fullfile(pwd,'TimeCells.mat'),'TodayTreadmillLog','TimeCells'); 
    outlines = cellfun(@bwboundaries,NeuronImage,'unif',0); 
    
    %Neuron centroids. 
    centroids = getNeuronCentroids(animal,date,session); 
        
    %Align and get indices where mouse was on treadmill. 
    [x,y,~,~,FToffset,~,~,time_interp] = AlignImagingToTracking(0.15,FT); 
    treadmillInds = getTreadmillEpochs(TodayTreadmillLog,time_interp); 
    nRuns = size(treadmillInds,1); 
  
%% Initialize movie properties.     
    %Tracking movie. 
    trackingread = dir('*.avi'); 
    trackingread = trackingread.name; 
       
    %Imaging movie.
    imagingmovie = VideoWriter(['TreadmillMovie_',movietype],'MPEG-4');
    imagingmovie.FrameRate = 20; 
    
    %Tracking movie.
    trackingread = VideoReader(trackingread); 
    trackingwrite = VideoWriter(['TrackingMovie_',movietype],'MPEG-4');
    trackingwrite.FrameRate = 20; 
    
    open(trackingwrite); 
    open(imagingmovie);
    
    tInc = 0;
     
    for thisEpoch=1:5
        if TodayTreadmillLog.complete(thisEpoch)
            sFrame = treadmillInds(thisEpoch,1) + FToffset;
            eFrame = treadmillInds(thisEpoch,2) + FToffset; 

            for i=sFrame:eFrame
%% Imaging movie. 
                %Get frame. 
                frame = h5read(h5file,'/Object',[1 1 i 1],[Xdim Ydim 1 1]);

                %Active neurons.
                active = find(FT(:,i));

                %Display the imaging movie frame. 
                imagesc(frame); caxis(clim); colormap gray;
                annotation(gcf,'textbox',[0.6, 0.8, 0.2, 0.07],'String',...
                {['t = ',num2str(round(tInc,1)),' seconds']},'Color','white','EdgeColor','white');
                annotation(gcf,'textbox',[0.2, 0.8, 0.1, 0.07],'String',...
                {['Lap ',num2str(thisEpoch)]},'Color','white','EdgeColor','white');
            
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
                clf;

%% Tracking movie
                %AVI file reads using time, not frames. Get the timestamp
                %from TodayTreadmillLog.
                trackingread.currentTime = TodayTreadmillLog.startts(thisEpoch) + tInc; 
                
                %Get frame. 
                frame = readFrame(trackingread); 
                
                %Display frame. 
                imagesc(frame); 
                annotation(gcf,'textbox',[0.6, 0.8, 0.2, 0.07],'String',...
                {['t = ',num2str(round(tInc,1)),' seconds']},'Color','red','EdgeColor','red');
                annotation(gcf,'textbox',[0.2, 0.8, 0.1, 0.07],'String',...
                {['Lap ',num2str(thisEpoch)]},'Color','red','EdgeColor','red');
                axis off; 
                
                %Write video to the tracking movie. 
                F = getframe(gcf);
                writeVideo(trackingwrite,F); 
                clf;
                
                %At the end of a lap, write a few blank frames. 
                if i==eFrame
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
end