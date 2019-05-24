function plotTimeCells(md,varargin)
%plotTimecells(MD,varargin)
%   
%   Plots single neuron responses in time during treadmill run. First
%   panel: spike-position trajectory. Second panel: Raster. Third panel:
%   Tuning curve. Press left/right arrows to scroll. Esc to exit. 
%
%   INPUTS
%       animal: Name of mouse (e.g., 'GCamp6f_45_treadmill').
%
%       date: Date of recording (e.g., '11_20_2015').
%
%       session: Session number. 
%   
%       T: Length of treadmill run. 
%

%% Grab inputs
    path = md.Location;
    cd(path);
    
    %Load time cell data. 
    try
        load(fullfile(pwd,'TimeCells.mat')); 
    catch
        [TimeCells,ratebylap,curves,movies,T,TodayTreadmillLog] = FindTimeCells(md,T); 
        tempInfo(md);
    end
    
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x)); 
    p.addParameter('dotplot',false,@(x) islogical(x)); 
    p.addParameter('placefield',false,@(x) islogical(x));
    p.addParameter('TimeCells',TimeCells,@(x) isnumeric(x)); 
    p.addParameter('singletraces',false,@(x) ischar(x));
    p.addParameter('plotit',true,@(x) islogical(x)); 
    
    p.parse(md,varargin{:});
    
    dotplot = p.Results.dotplot; 
    pf = p.Results.placefield; 
    TimeCells = p.Results.TimeCells;
    singletraces = p.Results.singletraces;
    plotit = p.Results.plotit;
    
%% 
    load('TemporalInfo.mat','Ispk');
    
    if pf            
        %Load place maps. 
        load(fullfile(path,'Placefields.mat'),'TMap_gauss','OccMap','pval'); 

        %Replace unoccupied bins with NaNs. 
        for i=1:length(TMap_gauss)
            TMap_gauss{i}(OccMap==0) = nan; 
        end
        
        pval = 1-pval; 
    end
    
    if singletraces
        try
            traces = load(fullfile(path,'TreadmillTraces.mat'),singletraces);
        catch
            TreadmillTraces(md);
            traces = load(fullfile(path,'TreadmillTraces.mat'),singletraces);
        end
    end
    
    %Extract the elements in structs. 
    if singletraces
        traces = traces.(singletraces);
    end
    x = movies.x;
    y = movies.y; 
    aviFrame = movies.t;
    FT = movies.FT; 
    delays = TodayTreadmillLog.delaysetting;
    complete = logical(TodayTreadmillLog.complete);
    alternation = strcmp(TodayTreadmillLog.direction,'alternation');

    %Basic setup. 
    [nNeurons,nFrames] = size(FT); 
    FT = logical(FT); 
    nBins = unique(sum(~isnan(ratebylap(delays==T & complete,:,1)),2));
    %pLaps = 0.2; 
    
    %Get indices for treadmill runs. 
    inds = getTreadmillEpochs(TodayTreadmillLog,aviFrame);
    temp = [];
    for thisLap=1:size(inds,1)
        if complete(thisLap) && delays(thisLap) == T
            temp = [temp,inds(thisLap,1):inds(thisLap,2)];
        end
    end
    treadmillruns = false(1,nFrames);
    treadmillruns(temp) = true; 

%% Plot setup. 
    thisNeuron = 1;
    keepgoing = 1;
    sf = 0.1;
    bins = [1:0.001:nBins]';
    t = linspace(0,T,length(bins));     %Time vector, smoothed.
    tCI = linspace(0,T,length(curves.ci{TimeCells(thisNeuron)}(1,:)));  %Time vector for CI interpolation.
    if singletraces
        tTraces = linspace(0,T,size(traces,2));
    end
    
    %Simplify the matrix and get rid of nans if any. 
    ratebylap = ratebylap(delays==T & complete,:,:);               %Laps that match delay duration. 
    ratebylap = ratebylap(:,~isnan(ratebylap(1,:,1)),:);
    if alternation
        turn = TodayTreadmillLog.choice(delays==T & complete); 
    end

%% For non-alternation.    
    if ~alternation
        while keepgoing
            %Smooth the tuning curve.  
            smoothfit = fit([1:nBins]',curves.tuning{TimeCells(thisNeuron)}','smoothingspline');
            curves.smoothed{TimeCells(thisNeuron)} = feval(smoothfit,bins);

            %Get confidence intervals interpolated.
            shuffmean = mean(curves.shuffle{TimeCells(thisNeuron)});
            CImean = interp1(tCI,shuffmean,t,'pchip');
            CIu = interp1(tCI,curves.ci{TimeCells(thisNeuron)}(1,:),t,'pchip');
            CIl = interp1(tCI,curves.ci{TimeCells(thisNeuron)}(2,:),t,'pchip');

            %Get time bins of significance.
            [SIGX,SIGY] = significance(t,curves.sig{TimeCells(thisNeuron)},...
                curves.smoothed{TimeCells(thisNeuron)},bins);

            %Plot. 
            if plotit
                f = figure(50);
                f.Position = [550 180 360 565];
                if dotplot
                    subplot(2,2,1);     %Dotplot. 
                        plot(x,y,x(treadmillruns & FT(TimeCells(thisNeuron),:)),y(treadmillruns & FT(TimeCells(thisNeuron),:)),'r.','MarkerSize',16);
                        axis off; title(['Neuron #',num2str(TimeCells(thisNeuron))]);
                    subplot(2,2,2);     %Raster. 
                        imagesc([0:T],[1:sum(complete)],ratebylap(:,:,TimeCells(thisNeuron)));
                            colormap gray; ylabel('Laps'); 
                elseif pf
                    subplot(2,2,1);     %Place map. 
                        h = imagesc(TMap_gauss{TimeCells(thisNeuron)});
                        set(h,'alphadata',~isnan(TMap_gauss{TimeCells(thisNeuron)}));
                            title(['p = ', num2str(pval(TimeCells(thisNeuron)))]);
                            axis off; colormap hot; freezeColors;
                    subplot(2,2,2);     %Raster. 
                        imagesc([0:T],[1:sum(complete)],ratebylap(:,:,TimeCells(thisNeuron)));
                            colormap gray; ylabel('Laps'); 
                            title(['Neuron #',num2str(TimeCells(thisNeuron))])
                else              
                    subplot(2,2,1:2);   %Raster. 
                        imagesc([0:T],[1:sum(complete)],ratebylap(:,:,TimeCells(thisNeuron)));
                            colormap gray; ylabel('Laps','fontsize',18); c = colorbar; c.Position(1) = 0.92;
                            title(['Neuron #',num2str(TimeCells(thisNeuron))]);
                            set(gca,'fontsize',16);
                end
                set(gca,'ytick',[1,sum(complete)]);

                subplot(2,2,3:4);   %Tuning curve.
                    yyaxis right; 
                    if singletraces
                        plot(tTraces,traces(:,:,TimeCells(thisNeuron)),...
                            '-','color',[.7 .7 .7 .2],'linewidth',2);
                            xlabel('Time [s]','fontsize',16); ylabel('\deltaF./\deltat');
                            set(gca,'fontsize',18,...
                                'ycolor',[.7 .7 .7],'linewidth',4,'tickdir','out');
                            ylim([min(min(traces(:,:,TimeCells(thisNeuron)))), ...
                                max(max(traces(:,:,TimeCells(thisNeuron))))]);
                    end

                    yyaxis left; 
                    hold on;
                    plot(t,curves.smoothed{TimeCells(thisNeuron)},'color',[0 .5 .5],...
                        'linewidth',5);
                    plot(t,CImean,'-b','linewidth',2);
                    plot(t,CIu,'--b',t,CIl,'--b');
                        Ylim = get(gca,'ylim');
                    plot(SIGX,SIGY+Ylim(2)*sf,'ro','linewidth',4);
                        xlim([0,T]);
                        ylabel('Rate','fontsize',16);
                        axis tight;
                        yLims = get(gca,'ylim');
                        ylim([0, yLims(2)*(1+sf)]);
                        set(gca,'ycolor',[0 .5 .5],'fontsize',15);
                    hold off;       
                    title(['I = ',num2str(round(Ispk(TimeCells(thisNeuron)),3)), ' bits']);

                %Scroll through neurons.
                [keepgoing,thisNeuron] = scroll(thisNeuron,length(TimeCells),f);
                close all;
            else 
                keepgoing = false;
            end
            
        end
%% Alternation
    else
        while keepgoing
            f = figure(50); 
            f.Position = [550 180 360 565];
                subplot(3,2,1:2);
                    plot(x,y,x(treadmillruns & FT(TimeCells(thisNeuron),:)),y(treadmillruns & FT(TimeCells(thisNeuron),:)),'r.','MarkerSize',16);
                    axis off; title(['Neuron #',num2str(TimeCells(thisNeuron))]);
            for lr=1:2
                good = turn==lr; 
                
                %Smooth the tuning curve.  
                smoothfit = fit([1:nBins]',curves.tuning{TimeCells(thisNeuron),lr}','smoothingspline');
                curves.smoothed{TimeCells(thisNeuron),lr} = feval(smoothfit,bins);

                %Get confidence intervals interpolated.
                shuffmean = mean(curves.shuffle{TimeCells(thisNeuron),lr});
                CImean = interp1(tCI,shuffmean,t,'pchip');
                CIu = interp1(tCI,curves.ci{TimeCells(thisNeuron),lr}(1,:),t,'pchip');
                CIl = interp1(tCI,curves.ci{TimeCells(thisNeuron),lr}(2,:),t,'pchip');
            
                %Get time bins of significance.
                [SIGX,SIGY] = significance(t,curves.sig{TimeCells(thisNeuron),lr},...
                    curves.smoothed{TimeCells(thisNeuron),lr},bins);
                           
                subplot(3,2,lr+2);
                    imagesc([0:T],[1:2:sum(delays(good)==T)],ratebylap(good,:,TimeCells(thisNeuron)));
                        colormap gray; ylabel('Laps');
                curveAX(lr) = subplot(3,2,lr+4);
                    plot(t,curves.smoothed{TimeCells(thisNeuron),lr},'-r','linewidth',2);
                    hold on; 
                    plot(t,CImean,'-b','linewidth',2);
                    plot(t,CIu,'--b',t,CIl,'--b');
                    Ylim = get(gca,'ylim');
                    xlim([0,T]);
                    
                    %Only plot significance asterisks if there was a
                    %response in one of the laps. 
                    if any(curves.smoothed{TimeCells(thisNeuron),lr})
                        plot(SIGX,SIGY+Ylim(2)*sf,'go');
                    end
                    
                    hold off;
                        xlabel('Time [s]'); ylabel('Rate'); 
                        yLims = get(gca,'ylim');
                        ylim([0, yLims(2)]);
                        set(gca,'ticklength',[0 0]);

                
            end
            
            %Normalize tuning curve axes. 
            curveYLims = [min([curveAX.YLim]), max([curveAX.YLim])];
            set(curveAX,'YLim',curveYLims); 
            
            %Scroll through neurons. 
            [keepgoing,thisNeuron] = scroll(thisNeuron,length(TimeCells),figure(50));
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