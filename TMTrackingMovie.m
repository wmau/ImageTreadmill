function TMTrackingMovie(md,varargin)
%
%
%

%%
    cd(md.Location);
    
    load(fullfile(pwd,'TimeCells.mat'),'TodayTreadmillLog','movies');
    load(fullfile(pwd,'FinalOutput.mat'),'FT'); 
    
    Pix2CM = md.Pix2CM; sf = 0.6246;
    [x,y,~,~,FToffset,~,~,~] = AlignImagingToTracking(Pix2CM,FT,10); 
    x = x./Pix2CM*sf; 
    y = y./Pix2CM*sf; 
    treadmillInds = getTreadmillEpochs(TodayTreadmillLog,movies.t); 
    nRuns = size(treadmillInds,1);
    
    %Tracking movie. 
    trackingread = dir('*.avi'); 
    trackingread = trackingread.name;
    trackingread = VideoReader(trackingread); 
    trackingwrite = VideoWriter(['TrackingMovie'],'MPEG-4');
    trackingwrite.FrameRate = 20; 
    
    open(trackingwrite);
    
    tfigure = figure('Position',[840 240 560 420]);
    
    tInc = 0; 
    for thisEpoch = 1:nRuns
        if TodayTreadmillLog.complete(thisEpoch)
            sFrame = treadmillInds(thisEpoch,1) + FToffset;
            eFrame = treadmillInds(thisEpoch,2) + FToffset; 
            
            for i=sFrame:eFrame
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
                
                if i==eFrame
                    clf; 
                    F = getframe(gcf); 
                    
                    for j=1:30
                        writeVideo(trackingwrite,F); 
                    end
                end
                
                tInc = tInc + 0.05;
            end
        end
        
        tInc = 0; 
    end
    
end