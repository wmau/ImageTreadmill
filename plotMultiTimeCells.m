function plotMultiTimeCells(batch_session_map,MD,Ts)
%
%
%

%% Organization. 
    %Partition the session data.
    nSessions = length(MD);
    dates = {MD.Date}; 
    paths = {MD.Location}; 
    animals = {MD.Animal};
    sessions = [MD.Session];
    
    %Preallocate.
    TIMECELLS = cell(nSessions,1);
    RATEBYLAP = cell(nSessions,1);
    CURVES = cell(nSessions,1);
    DELAYS = cell(nSessions,1); 
  
%% Gather all time cell data. 
    for i=1:nSessions
        cd(paths{i});
        
        try 
            load(fullfile(paths{i},'TimeCells.mat'),'TimeCells','ratebylap','curves','delays','T');
            
            if T~=Ts(i)
                error(['Delay duration specified in T(',num2str(i),') is different '...
                    'from the one saved for ',dates{i},'. Rerunning FindTimeCells '...
                    'using new T...']);
            end
            
            TIMECELLS{i} = TimeCells;
            RATEBYLAP{i} = ratebylap; 
            CURVES{i} = curves; 
            DELAYS{i} = delays; 
            
        catch
            [TIMECELLS{i},RATEBYLAP{i},CURVES{i},DELAYS{i}]...
                = FindTimeCells(animals{i},dates{i},sessions(i),Ts(i));
        end
        
    end
    
%% Find the indices in batch_session_map that correspond to the specified sessions. 
    regDates = {batch_session_map.session.date};
    regSessions = [batch_session_map.session.session];
      
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
    i = 1;
    sf = 0.1;
    pLaps = 0.25; 
    
    nBins = nan(nSessions,1);
    bins = cell(nSessions,1);
    tCI = cell(nSessions,1); 
    t = cell(nSessions,1); 
    critLaps = nan(nSessions,1);
    for thisSession=1:nSessions
        nBins(thisSession) = unique(sum(~isnan(RATEBYLAP{thisSession}(DELAYS{thisSession}==Ts(thisSession),:,1)),2));
        tCI{thisSession} = linspace(0,Ts(thisSession),length(CURVES{thisSession}.ci{1}(1,:)));        
        bins{thisSession} = [1:0.001:nBins(thisSession)]';
        t{thisSession} = linspace(0,Ts(thisSession),length(bins{thisSession}));
        
        nLaps = size(RATEBYLAP{thisSession}(DELAYS{thisSession}==Ts(thisSession),:,1),1);
        critLaps(thisSession) = round(nLaps*pLaps);
    end   
    
    f = figure('Position',[520    80   525   720]); 
    while keepgoing                   
        %Get the row index. 
        thisRow = uniqueRows(i);       
        neurons = MAP(thisRow,MAPinds);     %Neurons in this row. 

        %For each session, plot its ratebylap. If subplot has a title,
        %neurno was mapped, but inactive on stem. 
        for thisSession=1:nSessions
            n = neurons(thisSession);
            
            %Raster. 
            rasterAX(thisSession) = subplot(nSessions,2,thisSession*2-1);
            if n  
                imagesc([0:Ts(thisSession)],[1:5:sum(DELAYS{thisSession}==Ts(thisSession))],...
                    RATEBYLAP{thisSession}(DELAYS{thisSession}==Ts(thisSession),:,n));
                colormap gray; ylabel('Laps'); title(['Neuron #',num2str(n)]);
            else
                imagesc(0); axis off; 
            end
            
            %Tuning curve. 
            curveAX(thisSession) = subplot(nSessions,2,thisSession*2);
            if n
                %Smooth tuning curve.                
                smoothfit = fit([1:nBins(thisSession)]',CURVES{thisSession}.tuning{n}','smoothingspline');
                CURVES{thisSession}.smoothed{n} = feval(smoothfit,bins{thisSession});
                
                %Confidence interval interpolation.               
                shuffmean = mean(CURVES{thisSession}.shuffle{n});
                CImean = interp1(tCI{thisSession},shuffmean,t{thisSession},'pchip');
                CIu = interp1(tCI{thisSession},CURVES{thisSession}.ci{n}(1,:),t{thisSession},'phcip');
                CIl = interp1(tCI{thisSession},CURVES{thisSession}.ci{n}(2,:),t{thisSession},'phcip');
                
                [SIGX,SIGY] = significance(t{thisSession},CURVES{thisSession}.sig{n},...
                    CURVES{thisSession}.smoothed{n},bins{thisSession});
                
                %Plot tuning curve and confidence intervals. 
                plot(t{thisSession},CURVES{thisSession}.smoothed{n},'-r','linewidth',2);
                hold on;
                plot(t{thisSession},CImean,'-b','linewidth',2);
                plot(t{thisSession},CIu,'--b',t{thisSession},CIl,'--b');  
                Ylim = get(gca,'ylim');
        
                %If there are enough laps, plot significance asterisks. 
                if sum(any(RATEBYLAP{thisSession}(:,:,n),2)) > critLaps(thisSession)
                    plot(SIGX,SIGY+Ylim(2)*sf,'r*');                   
                end
                
                %Labels. 
                xlabel('Time [s]'); ylabel('Rate');
                yLims = get(gca,'ylim');
                ylim([0,yLims(2)]); xlim([0,t{thisSession}(end)]);
                set(gca,'ticklength',[0 0]);
                hold off; 
            else
                plot(t{thisSession},zeros(length(t{thisSession}),1),'-r','linewidth',2); 
                yLims = get(gca,'ylim'); xlim([0,t{thisSession}(end)]);
                ylim([0,yLims(2)]);
            end
        end
        
        rasterXLims = [min([rasterAX.XLim]), max([rasterAX.XLim])];
        curveXLims = [min([curveAX.XLim]), max([curveAX.XLim])];
        curveYLims = [min([curveAX.YLim]), max([curveAX.YLim])];
        set(rasterAX,'XLim',rasterXLims);
        set(curveAX,'XLim',curveXLims,'YLim',curveYLims);
        
        figure(f);
            [~,~,key] = ginput(1);
            if key == 29 && i < length(uniqueRows)
                i = i + 1; 
            elseif key == 28 && i ~= 1
                i = i - 1; 
            elseif key == 27
                keepgoing = 0; 
                close(figure(f)); 
            end

    end
    
end
        
function [SIGX,SIGY] = significance(t,sigcurve,smoothedcurve,bins)
    %Get indices of significance.
    sigT = find(sigcurve);
    
    %Find the corresponding index in the smoothed bins.
    [~,inds] = ismember(sigT,bins);
    
    %Values to plot.
    SIGX = t(inds);
    SIGY = smoothedcurve(inds);
end