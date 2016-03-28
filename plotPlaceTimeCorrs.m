function plotPlaceTimeCorrs(MAPMD,MD1,MD2,Ts)
%
%
%

%%
    [pCorr,tCorr,MAP,MAPcols,DATA,noi] = PlaceTimeCorr(MAPMD,MD1,MD2); 
    RATEBYLAP = DATA.ratebylap; 
    CURVES = DATA.curves;
    DELAYS = DATA.delays; 
    COMPLETE = DATA.complete; 
    PFS = DATA.placefields; 
    PVALS = DATA.placefieldpvals;
    OCCMAPS = DATA.occmaps;
    
    MD = [MD1 MD2];
    dates = {MD.Date};
    dateTitles = dates; 
    for s=1:2
        dateTitles{s}(3:3:6) = '-';
    end
    
    keepgoing = 1;
    sf = 0.1; 
    pLaps = 0.25;
    rows = find(ismember(MAP(:,MAPcols(1)),noi));
    i = 1;
    
    %Get the basics for each session like confidence interval allocations
    %and delay duration. 
    nBins = nan(2,1);       %Number of bins corresponding to delay duration.
    bins = cell(2,1);       %Linear space for smoothing.
    tCI = cell(2,1);        %Linear space for interpolating CIs. 
    t = cell(2,1);          %Time vector for the x-axis of smoothed tuning curve. 
    pfExist = logical(zeros(2,1)); 
    critLaps = nan(2,1);
    
    for s=1:2
        nBins(s) = max(sum(~isnan(RATEBYLAP{s}(DELAYS{s}==Ts(s),:,1)),2));
        tCI{s} = linspace(0,Ts(s),length(CURVES{s}.ci{1}(1,:)));        
        bins{s} = [1:0.001:nBins(s)]';
        t{s} = linspace(0,Ts(s),length(bins{s}));
        
        nLaps = size(RATEBYLAP{s}(DELAYS{s}==Ts(s),:,1),1);
        critLaps(s) = round(nLaps*pLaps);         %Critical number of laps. 
    end   
    
    %Main plotting. 
    f = figure('Position',[-1300, -40, 520, 360]); 
    while keepgoing                   
        r = rows(i);
        neurons = MAP(r,MAPcols);
        cmax = zeros(1,2); 
        
        for s=1:2
            n = neurons(s); 
            
            %PLACE FIELDS. 
            if logical(all(isnan(PFS{s}{n}(:))))
                pfAX(s) = subplot(2,3,s*3-2);
                    h = imagesc(OCCMAPS{s}); 
                    set(h,'alphadata',~isnan(OCCMAPS{s}));
                    axis off; colormap gray; freezeColors;
                    pfExist(s) = false; 
            else
                pfAX(s) = subplot(2,3,s*3-2);
                    h = imagesc(PFS{s}{n}); 
                    set(h,'alphadata',~isnan(PFS{s}{n}));
                    axis off; colormap hot; freezeColors;

                cmax(s) = max(PFS{s}{n}(:));

                if s==1
                    title({['p=',num2str(PVALS{s}(n))],...
                       ['Corr p=',num2str(pCorr(n,2))]});
                else
                    title(['p=',num2str(PVALS{s}(n))]);
                end
                
                pfExist(s) = true; 
            end
            
            %RASTER. 
            rasterAX(s) = subplot(2,3,s*3-1); 
            runningAtT = DELAYS{s}==Ts(s);
            goodLaps = runningAtT & COMPLETE{s}; 
            
            plotme = RATEBYLAP{s}(DELAYS{s}==Ts(s) & COMPLETE{s},:,n);
            plotme = plotme(:,~isnan(plotme(1,:)));
            
            imagesc(0:Ts(s),1:5:sum(goodLaps),plotme);
            colormap gray; freezeColors;
            ylabel('Laps'); title(['Neuron #',num2str(n)]); 
            
            %TUNING CURVE
            curveAX(s) = subplot(2,3,s*3); 
            
            smoothfit = fit([1:nBins(s)]',CURVES{s}.tuning{n}','smoothingspline');
            CURVES{s}.smoothed{n} = feval(smoothfit,bins{s});

            %Confidence interval interpolation.               
            shuffmean = mean(CURVES{s}.shuffle{n});
            CImean = interp1(tCI{s},shuffmean,t{s},'pchip');
            CIu = interp1(tCI{s},CURVES{s}.ci{n}(1,:),t{s},'phcip');
            CIl = interp1(tCI{s},CURVES{s}.ci{n}(2,:),t{s},'phcip');

            %Plot tuning curve and confidence intervals. 
            plot(t{s},CURVES{s}.smoothed{n},'-r','linewidth',2);
            hold on;
            plot(t{s},CImean,'-b','linewidth',2);
            plot(t{s},CIu,'--b',t{s},CIl,'--b');  
            Ylim = get(gca,'ylim');

            %If there are enough laps, plot significance asterisks. 
            if sum(any(RATEBYLAP{s}(:,:,n),2)) > critLaps(s)
                %Significance asterisks. 
                [SIGX,SIGY] = significance_asterisks(t{s},CURVES{s}.sig{n},...
                    CURVES{s}.smoothed{n},bins{s});

                plot(SIGX,SIGY+Ylim(2)*sf,'r*');                   
            end
            
            %Labels. 
            if s==1
                title({dateTitles{s},...
                    ['Corr p=',num2str(tCorr(n,2))]});
            else
                title(dateTitles{s});
            end
            
            xlabel('Time [s]'); ylabel('Rate');
            yLims = get(gca,'ylim');
            ylim([0,yLims(2)]); xlim([0,t{s}(end)]);
            set(gca,'ticklength',[0 0]);
            hold off; freezeColors;
            
           
        end
        
        %Normalize the axes. 
        rasterXLims = [min([rasterAX.XLim]), max([rasterAX.XLim])];
        curveXLims = [min([curveAX.XLim]), max([curveAX.XLim])];
        curveYLims = [min([curveAX.YLim]), max([curveAX.YLim])];

        %Normalizing axes. 
        for ss=1:2
            if pfExist(ss)
                subplot(2,3,ss*3-2);
                unfreezeColors; colormap hot; 
            end
        end
        clims = [pfAX.CLim];
        pfLims = [0, max(clims(clims~=1))];
        set(pfAX(pfExist),'CLim',pfLims); 

        for ss=1:2
            subplot(2,3,ss*3-2);
            freezeColors;
        end

        set(rasterAX,'XLim',rasterXLims);
        set(curveAX,'XLim',curveXLims,'YLim',curveYLims)
        
        %Scroll through neurons. 
        [keepgoing,i] = scroll(i,length(rows),f);
    end
    
end