function plotMultiTimeCells(mapMD,MD,Ts)
%plotMultiTimeCells(MAPlocation,MD,Ts)
%
%   Plots time cells across multiple sessions. Use left and right arrow
%   keys to scroll, esc to exit. 
%
%   INPUTS
%       MAPlocation: Directory containing struct from neuron_batch_reg().
%
%       MD: MasterDirectory entries. 
%
%       Ts: Vector the same size as MD specifying the delay durations for
%       each session. 
%

%% Organization. 
    initDir = pwd; 
    load(fullfile(mapMD.Location),'batch_session_map.mat'); 

    %Partition the session data.
    nSessions = length(MD);
    dates = {MD.Date}; 
    sessions = [MD.Session];
    
    %Setting up for titling the plots. 
    dateTitles = dates;
    for i=1:nSessions
        dateTitles{i}(3:3:6) = '-';
    end
  
%% Gather all time cell data. 
    [TIMECELLS,RATEBYLAP,CURVES,DELAYS,COMPLETE] = CompileTimeCellData(MD,Ts);
    
%% Find the indices in batch_session_map that correspond to the specified sessions. 
    regDates = {batch_session_map.session.Date};
    regSessions = [batch_session_map.session.Session];
      
    %Eliminate first column to match indices.
    MAP = batch_session_map.map(:,2:end);           
    
    %Get column indices of batch_session_map.map. 
    MAPinds = zeros(nSessions,1); 
    rows = cell(nSessions,1); 
    for i=1:nSessions
        try
            MAPinds(i) = find(ismember(regDates,dates{i}) & ismember(regSessions,sessions(i)));
        catch
            error(['Error in above. Possible reason: MD input has not been registered yet. '...
                'Run neuron_reg_batch...']);
        end
        
        %Get row indices of MAP containing time cells in at least one
        %session. 
        rows{i} = find(ismember(MAP(:,MAPinds(i)),TIMECELLS{i}));
    end
    
    %All the row indices of MAP that have time cells. 
    uniqueRows = unique(cell2mat(rows)); 
    
%% Plot. 
    keepgoing = 1;                  
    i = 1;                          %Row of neuron map. 
    sf = 0.1;                       %For significance asterisks. 
    pLaps = 0.20;                   %Proportion of laps needed to be significant. 
    
    %Get the basics for each session like confidence interval allocations
    %and delay duration. 
    nBins = nan(nSessions,1);       %Number of bins corresponding to delay duration.
    bins = cell(nSessions,1);       %Linear space for smoothing.
    tCI = cell(nSessions,1);        %Linear space for interpolating CIs. 
    t = cell(nSessions,1);          %Time vector for the x-axis of smoothed tuning curve. 
    critLaps = nan(nSessions,1);
    for thisSession=1:nSessions
        nBins(thisSession) = max(sum(~isnan(RATEBYLAP{thisSession}(DELAYS{thisSession}==Ts(thisSession),:,1)),2));
        tCI{thisSession} = linspace(0,Ts(thisSession),length(CURVES{thisSession}.ci{1}(1,:)));        
        bins{thisSession} = [1:0.001:nBins(thisSession)]';
        t{thisSession} = linspace(0,Ts(thisSession),length(bins{thisSession}));
        
        nLaps = size(RATEBYLAP{thisSession}(DELAYS{thisSession}==Ts(thisSession),:,1),1);
        critLaps(thisSession) = round(nLaps*pLaps);         %Critical number of laps. 
    end   
    
    %Main plotting. 
    f = figure('Position',[520    80   525   720]); 
    while keepgoing                   
        %Get the row index. 
        thisRow = uniqueRows(i);       
        neurons = MAP(thisRow,MAPinds);     %Neurons in this row. 

        %For each session, plot its ratebylap. If subplot has a title,
        %neuron was mapped, but inactive on stem. 
        for thisSession=1:nSessions
            n = neurons(thisSession);
            
            %Raster. 
            rasterAX(thisSession) = subplot(nSessions,2,thisSession*2-1);
            if ~isnan(n) && n~=0
                runningAtT = DELAYS{thisSession}==Ts(thisSession);      %Logical. Is the mouse running for Ts this lap?
                goodLaps = runningAtT & COMPLETE{thisSession};          %Logical. Is the mouse running for Ts for the whole time this lap?
                
                %Restrict plotting to complete laps with the mouse running
                %at the specified T.              
                plotme = RATEBYLAP{thisSession}(DELAYS{thisSession}==Ts(thisSession) & COMPLETE{thisSession},:,n);
                plotme = plotme(:,~isnan(plotme(1,:)));

                %Plot raster. 
                imagesc(0:Ts(thisSession),1:5:sum(goodLaps),plotme);
                colormap gray; ylabel('Laps'); title(['Neuron #',num2str(n)]);
            else
                imagesc(0); axis off; 
            end
            
            %Tuning curve. 
            curveAX(thisSession) = subplot(nSessions,2,thisSession*2);
            if ~isnan(n) && n~=0
                %Smooth tuning curve.                
                smoothfit = fit([1:nBins(thisSession)]',CURVES{thisSession}.tuning{n}','smoothingspline');
                CURVES{thisSession}.smoothed{n} = feval(smoothfit,bins{thisSession});
                
                %Confidence interval interpolation.               
                shuffmean = mean(CURVES{thisSession}.shuffle{n});
                CImean = interp1(tCI{thisSession},shuffmean,t{thisSession},'pchip');
                CIu = interp1(tCI{thisSession},CURVES{thisSession}.ci{n}(1,:),t{thisSession},'phcip');
                CIl = interp1(tCI{thisSession},CURVES{thisSession}.ci{n}(2,:),t{thisSession},'phcip');
                               
                %Plot tuning curve and confidence intervals. 
                plot(t{thisSession},CURVES{thisSession}.smoothed{n},'-r','linewidth',2);
                hold on;
                plot(t{thisSession},CImean,'-b','linewidth',2);
                plot(t{thisSession},CIu,'--b',t{thisSession},CIl,'--b');  
                Ylim = get(gca,'ylim');
        
                %If there are enough laps, plot significance asterisks. 
                if sum(any(RATEBYLAP{thisSession}(:,:,n),2)) > critLaps(thisSession)
                    %Significance asterisks. 
                    [SIGX,SIGY] = significance_asterisks(t{thisSession},CURVES{thisSession}.sig{n},...
                        CURVES{thisSession}.smoothed{n},bins{thisSession});
                    
                    plot(SIGX,SIGY+Ylim(2)*sf,'r*');                   
                end
                
                %Labels. 
                title(dateTitles{thisSession});
                xlabel('Time [s]'); ylabel('Rate');
                yLims = get(gca,'ylim');
                ylim([0,yLims(2)]); xlim([0,t{thisSession}(end)]);
                set(gca,'ticklength',[0 0]);
                hold off; 
            else
                %Flat line. 
                plot(t{thisSession},zeros(length(t{thisSession}),1),'-r','linewidth',2);
                title(dateTitles{thisSession});
                yLims = get(gca,'ylim'); xlim([0,t{thisSession}(end)]);
                ylim([0,yLims(2)]);
            end
        end
        
        %Normalize the axes. 
        rasterXLims = [min([rasterAX.XLim]), max([rasterAX.XLim])];
        curveXLims = [min([curveAX.XLim]), max([curveAX.XLim])];
        curveYLims = [min([curveAX.YLim]), max([curveAX.YLim])];
        set(rasterAX,'XLim',rasterXLims);
        set(curveAX,'XLim',curveXLims,'YLim',curveYLims);
        
        %Scroll through neurons. 
        [keepgoing,i] = scroll(i,length(uniqueRows),f);

    end
    
    cd(initDir); 
end