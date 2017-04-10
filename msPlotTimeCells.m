function msPlotTimeCells(md,varargin)
%msPlotTimeCells(MD,Ts,varargin)
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

%% Parse inputs.
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x));
    p.addParameter('timecells',[]);
    p.addParameter('pf',false,@(x) islogical(x));
    
    p.parse(md,varargin{:});
    timecells = p.Results.timecells;
    pf = p.Results.pf;

%% Organization. 
    initDir = pwd; 
    mapMD = getMapMD(md(1));
    load(fullfile(mapMD.Location,'batch_session_map.mat')); 
     
    if pf, nCols=3;  
    else nCols=2; end
    
    %Partition the session data.
    nSessions = length(md);
    dates = {md.Date}; 
    sessions = [md.Session];
    
    %Setting up for titling the plots. 
    dateTitles = dates;
    for i=1:nSessions
        dateTitles{i}(3:3:6) = '-';
    end
  
%% Gather all time cell data. 
    args = {'timecells',...
            't',...
            'ratebylap',...
            'curves',...
            'delays',...
            'complete',...
            'ti',...
            'tfcorr'};
    if pf
        args{end+1} = 'placefields'; 
        args{end+1} = 'runoccmaps'; 
        args{end+1} = 'placefieldpvals';
        args{end+1} = 'ttl';
    end
        
    DATA = CompileMultiSessionData(md,args);
    
    TIMECELLS = DATA.timecells; 
    RATEBYLAP = DATA.ratebylap; 
    CURVES = DATA.curves; 
    DELAYS = DATA.delays; 
    COMPLETE = DATA.complete; 
    TI = DATA.ti;
    TFCORR = DATA.tfcorr;
    Ts = cell2mat(DATA.t); 
    if pf
        PFS = DATA.placefields; 
        RUNOCCMAPS = DATA.runoccmaps; 
        PVALS = DATA.placefieldpvals; 
        TTL = DATA.ttl;
    end
    
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
            error(['Error in above. Possible reason: MD input ', num2str(i),...
                ' has not been registered yet. Run neuron_reg_batch...']);
        end
        
        %Get row indices of MAP containing time cells in at least one
        %session. 
        rows{i} = find(ismember(MAP(:,MAPinds(i)),TIMECELLS{i}));
    end
 
    if ~isempty(timecells)
        %Only look at select cells.
        rows = find(ismember(MAP(:,MAPinds(1)),timecells));
    else 
        %All the row indices of MAP that have time cells. 
        rows = unique(cell2mat(rows));
    end
    
%% Plot. 
    keepgoing = 1;                  
    i = 1;                          %Row of neuron map. 
    sf = 0.1;                       %For significance asterisks. 
    pLaps = 0.25;                   %Proportion of laps needed to be significant. 
    
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
    if pf 
        fPos = [-1300 -40 630 180*nSessions];
    else 
        fPos = [-1300 -40 370 180*nSessions];
    end
    
    
    while keepgoing        
        f = figure('Position',fPos); 
        
        %Get the row index. 
        thisRow = rows(i);       
        neurons = MAP(thisRow,MAPinds);     %Neurons in this row. 
        cmax = zeros(1,nSessions); 
        pfExist = false(1,nSessions); 

        %For each session, plot its ratebylap. If subplot has a title,
        %neuron was mapped, but inactive on stem. 
        for thisSession=1:nSessions
            n = neurons(thisSession);
            
            if isnan(n) || n==0, detection = ' Not detected'; 
            else, detection = []; end
            
            %PLACE FIELD. 
            if ~isnan(n) && n~=0 && pf
                cd(md(thisSession).Location); 
                load('Pos_align.mat','x_adj_cm','y_adj_cm','xmin','ymin','xmax','ymax');
                x = ymin:ymax; 
                y = xmin:xmax; 
                
                %Get sections. 
                direction = TTL{thisSession}.direction;         
                sections = sections_treadmill(x_adj_cm,y_adj_cm,direction,false); 
                
                f.Position = fPos;
                	pfAX(thisSession) = subplot(nSessions,nCols,thisSession*nCols-2);
                    h = imagesc(x,y,PFS{thisSession}{n}); 
                    set(h,'alphadata',~isnan(PFS{thisSession}{n}));
                    axis off; colormap hot; freezeColors; 
                    hold on; 
                    
                    rectangle(  'position',[sections.center.y(1),...
                                sections.center.x(1),...
                                sections.center.y(3)-sections.center.y(2),...
                                sections.center.x(2)-sections.center.x(1)],...
                                'edgecolor',[139 69 19]./255,...
                                'linewidth',2,...
                                'linestyle','--');

                    
                    %Peak place field. 
                    cmax(thisSession) = max(PFS{thisSession}{n}(:));            
                    
                    %P-value of the place field. 
                    title(['p=',num2str(PVALS{thisSession}(n))]);
                    
                    pfExist(thisSession) = true; 
                    
                    colorbar;
            end
            
            %Occupancy map. 
            if ~isnan(n) && n~=0 && pf && logical(all(PFS{thisSession}{n}(:)==0))
                f.Position = [-1300 -40 520 180*nSessions];
                    pfAX(thisSession) = subplot(nSessions,nCols,thisSession*nCols-2);
                    h = imagesc(RUNOCCMAPS{thisSession}); 
                    set(h,'alphadata',~isnan(RUNOCCMAPS{thisSession}));
                    axis off; colormap gray; freezeColors;
                    
                pfExist(thisSession) = false; 
            end
            
            if (isnan(n) || n==0) && pf     %Blank.
                f.Position = fPos;
                    pfAX(thisSession) = subplot(nSessions,nCols,thisSession*nCols-2);
                    imagesc(0); axis off; colormap gray; freezeColors;
                
                pfExist(thisSession) = false; 
            end
            
            %RASTER. 
            rasterAX(thisSession) = subplot(nSessions,nCols,thisSession*nCols-1);
            if ~isnan(n) && n~=0
                runningAtT = DELAYS{thisSession}==Ts(thisSession);      %Logical. Is the mouse running for Ts this lap?
                goodLaps = runningAtT & COMPLETE{thisSession};          %Logical. Is the mouse running for Ts for the whole time this lap?
                
                %Restrict plotting to complete laps with the mouse running
                %at the specified T.              
                plotme = RATEBYLAP{thisSession}(DELAYS{thisSession}==Ts(thisSession) & COMPLETE{thisSession},:,n);
                plotme = plotme(:,~isnan(plotme(1,:)));

                %Plot raster. 
                imagesc(0:Ts(thisSession),1:sum(goodLaps),plotme);
                colormap gray; freezeColors;
                
                if thisSession~=nSessions, set(gca,'xtick',[]); end
                
                ylabel('Laps'); title(['Neuron #',num2str(n)]); 
            else
                imagesc(0); 
                colormap gray; axis off; freezeColors
            end
            set(gca,'ytick',[1,sum(COMPLETE{thisSession})]);
            
            %TUNING CURVE. 
            curveAX(thisSession) = subplot(nSessions,nCols,thisSession*nCols);
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
                plot(t{thisSession},CURVES{thisSession}.smoothed{n},'color',[0 .5 .5],...
                    'linewidth',2);
                hold on;
                plot(t{thisSession},CImean,'-b','linewidth',2);
                plot(t{thisSession},CIu,'--b',t{thisSession},CIl,'--b');  
                Ylim = get(gca,'ylim');
        
                %If there are enough laps, plot significance asterisks. 
                if ismember(n,TIMECELLS{thisSession})
                    %Significance asterisks. 
                    [SIGX,SIGY] = significance_asterisks(t{thisSession},CURVES{thisSession}.sig{n},...
                        CURVES{thisSession}.smoothed{n},bins{thisSession});
                    
                    plot(SIGX,SIGY+Ylim(2)*sf,'ro','linewidth',2);                   
                end
                
                %Labels. 
                title(dateTitles{thisSession});
                yLims = get(gca,'ylim');
                ylim([0,yLims(2)]); xlim([0,t{thisSession}(end)]);
                ylabel('Rate');
                
                if thisSession~=nSessions, set(gca,'xtick',[]); end
                if thisSession==nSessions, xlabel('Time [s]'); end
                
                set(gca,'tickdir','out','linewidth',4);
                hold off; freezeColors
            else
                %Flat line. 
                plot(t{thisSession},zeros(length(t{thisSession}),1),'-r','linewidth',2);
                title([dateTitles{thisSession} detection]);
                set(gca,'xtick',[]);
                yLims = get(gca,'ylim'); xlim([0,t{thisSession}(end)]);
                ylim([0,yLims(2)]); freezeColors
            end
        end
        
        %Normalize the axes. 
        rasterXLims = [min([rasterAX.XLim]), max([rasterAX.XLim])];
        curveXLims = [min([curveAX.XLim]), max([curveAX.XLim])];
        curveYLims = [min([curveAX.YLim]), max([curveAX.YLim])];
        
        %Label temporal information.
        for thisSession=1:nSessions
            n = neurons(thisSession);
            
            if ~isnan(n) && n~=0
                ax = subplot(nSessions,nCols,thisSession*nCols);
                set(ax,'units','normalized');
                text(2.5,0.9*curveYLims(2),['TI = ',num2str(round(TI{thisSession}(n),3)), ' bits']); 
                
                if thisSession~=nSessions
                    TCs = getTimeCells(md(thisSession));
                    TCs = EliminateUncertainMatches([md(thisSession),md(thisSession+1)],TCs);
                    
                    crit = 0.01/length(TCs);
                    if TFCORR{thisSession}(n,2) < crit
                        c = [0 .5 .5]; 
                    else, c = 'r';
                    end 
                    text(6,-.5,['R = ',num2str(round(TFCORR{thisSession}(n,1),3))],'color',c)     
                end
            end
     
        end
        
        if pf 
            for j=1:nSessions
                if pfExist(j)
                    subplot(nSessions,nCols,j*nCols-2);
                    unfreezeColors; colormap hot; 
                end
            end
            clims = [pfAX.CLim];
            pfLims = [0, max(clims(clims~=1))];
            try set(pfAX(pfExist),'CLim',pfLims); 
            catch, set(pfAX(pfExist),'CLim',[0 1]); end
            for j=1:nSessions
                subplot(nSessions,nCols,j*nCols-2);
                freezeColors;
            end
        end  
        
        set(rasterAX,'XLim',rasterXLims);
        set(curveAX,'XLim',curveXLims,'YLim',curveYLims);
        
        %Scroll through neurons. 
        [keepgoing,i] = scroll(i,length(rows),f);
        close all;
    end
    
    cd(initDir); 
end