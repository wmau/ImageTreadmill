function [raster,smoothCurve,curve] = LinearizedPF_raster(md,varargin)
%
%
%

%% Parse inputs.
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x));
    p.addParameter('plotit',true,@(x) islogical(x));
    p.addParameter('minspeed',3,@(x) isscalar(x));
    p.addParameter('nBins',80,@(x) isscalar(x)); 
    p.addParameter('neurons',[],@(x) isnumeric(x));
    p.addParameter('saveBool',false,@(x) islogical(x));
    
    p.parse(md,varargin{:});
    
    plotit = p.Results.plotit;
    minspeed = p.Results.minspeed; 
    nBins = p.Results.nBins; 
    neurons = p.Results.neurons; 
    saveBool = p.Results.saveBool;
    
    if isempty(neurons)
        neurons = getPlaceCells(md,.01);
    end
    nNeurons = length(neurons);
    
%%
    currDir = pwd;
    cd(md.Location); 
    
    %Get treadmill log for excluding treadmill epochs. 
    load('TimeCells.mat','TodayTreadmillLog'); 
    d = TodayTreadmillLog.direction; 
    
    %Find direction for linearizing trajectory. 
    if strfind(d,'left')
        mazetype = 'left';
    elseif strfind(d,'right')
        mazetype = 'right';
    elseif strfind(d,'alternation')
        mazetype = 'tmaze';
    end
    
       %Load aligned position data. 
    load(fullfile(pwd,'Pos_align.mat'),'x_adj_cm','y_adj_cm','speed','time_interp','PSAbool');
    x=x_adj_cm; y=y_adj_cm; PSAbool=logical(PSAbool); clear x_adj_cm y_adj_cm;
    [~,nFrames] = size(PSAbool); 
    
    %Exclude treadmill epochs. 
    inds = TodayTreadmillLog.inds;
    excludeFrames=[]; 
    nSeconds = 2;
    extraExclude = 20*nSeconds;
    for e=1:size(inds,1)
        window = (inds(e,1)-extraExclude):(inds(e,2)+extraExclude);
        window = window(window>1);
        window = window(window<nFrames);

        excludeFrames = [excludeFrames, window];
    end
    onTM = false(1,nFrames); 
    onTM(excludeFrames) = true; 

    %Speed threshold. 
    isrunning = speed>minspeed; 
    
%% Linearize trajectory and bin responses spatially.
    %Linearized trajectory. 
    X = LinearizeTrajectory_treadmill(x,y,mazetype); 
    parsed = postrials_treadmill(md);
    trials = parsed.summary(:,1)'; 
    
    keepgoing = true;
    bins = [1:0.001:nBins]';
    raster = zeros(max(trials),nBins,nNeurons);
    curve = zeros(nNeurons,nBins);
    smoothCurve = zeros(nNeurons,size(bins,1)); 
    
    for thisNeuron = 1:nNeurons    
        for thisTrial = trials
            if ~isnan(neurons(thisNeuron)) && neurons(thisNeuron)>0
                spkpos = X(parsed.trial==thisTrial & PSAbool(neurons(thisNeuron),:) & isrunning & ~onTM);
            else
                spkpos = [];
            end

            [occ,edges] = histcounts(X(parsed.trial==thisTrial),nBins);
            binned = histcounts(spkpos,edges);

            raster(thisTrial,:,thisNeuron) = binned./occ;

        end
        
        %Smooth.
        curve(thisNeuron,:) = nanmean(raster(:,:,thisNeuron));
        smoothfit = fit([1:nBins]',curve(thisNeuron,:)','smoothingspline');
        smoothCurve(thisNeuron,:) = feval(smoothfit,bins)';
    end
        
    %Plot.
    if plotit
        thisNeuron = 1;
        while keepgoing

        f = figure(49); 
        f.Position = [550 180 360 565];
        subplot(2,1,1); 
            imagesc(raster(:,:,thisNeuron)); 
            colormap hot; caxis([0 1.2]);
            set(gca,'xtick',[],'ytick',[1,max(trials)],'linewidth',4,...
                'fontsize',12);
            ylabel('Laps','fontsize',15); 
        subplot(2,1,2); 
            plot(bins,smoothCurve(thisNeuron,:),'color',[.58 .44 .86],'linewidth',5);
            xlabel('Linearized Distance','fontsize',15);
            ylabel('Rate','fontsize',15);
            yLims = get(gca,'ylim');
            ylim([0 yLims(2)]);
            set(gca,'linewidth',4,'tickdir','out','fontsize',12);

        [keepgoing,thisNeuron] = scroll(thisNeuron,nNeurons,f);
        end
    end
    
    cd(md.Location);
    if saveBool
        save('SpatialTraces.mat','raster','curve','smoothCurve');
    end
    
    cd(currDir);
    
end