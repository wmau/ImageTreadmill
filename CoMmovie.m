function CoMmovie(animal,date,session)
%
%
%

%% 
    close all;
    ChangeDirectory(animal,date,session);
    
    load(fullfile(pwd,'ProcOut_minlength_2.mat'),'FT','NeuronImage','Xdim','Ydim');
    load(fullfile(pwd,'TimeCells.mat'),'TodayTreadmillLog');
    
    centroids = getNeuronCentroids(animal,date,session,2); 
    
    [x,y,~,~,FToffset,~,~,time_interp] = AlignImagingToTracking(0.15,FT); 
    treadmillInds = getTreadmillEpochs(TodayTreadmillLog,time_interp); 
    nRuns = size(treadmillInds,1); 
    
    CoMmovie = VideoWriter('CoMmovie','MPEG-4'); 
    CoMmovie.FrameRate = 20; 
    
    open(CoMmovie); 
    
    %Set up colors.
    colors = rand(nRuns,3); 
    
    figure; hold on;
    for thisEpoch=1:nRuns
        if TodayTreadmillLog.complete(thisEpoch)
            sFrame = treadmillInds(thisEpoch,1) + FToffset;
            eFrame = treadmillInds(thisEpoch,2) + FToffset; 
            
            for i=sFrame:eFrame
                active = find(FT(:,i)); 

                FTCoMx = mean(centroids(active,1)); 
                FTCoMy = mean(centroids(active,2)); 
                
                if i==sFrame
                    plot(FTCoMx,FTCoMy,'.','color',colors(thisEpoch,:));
                else
                    line([lastx FTCoMx],[lasty FTCoMy],'color',colors(thisEpoch,:),'linewidth',5);
                end
                
                lastx = FTCoMx; lasty = FTCoMy; 
                
                axis equal; axis off; 
                F = getframe(gcf); 
                writeVideo(CoMmovie,F);
                clear F;
            end
        end
    end
    hold off; 
    
end