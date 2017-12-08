function [raster,smoothCurve,curve,X,binocc,parsed] = ...
    LinearizedPF_raster(md,varargin)
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
    p.addParameter('shuffle',true,@(x) islogical(x)); 
    p.addParameter('noTreadmill',true,@(x) islogical(x)); 
    p.addParameter('plotTrials',true,@(x) islogical(x));
    
    p.parse(md,varargin{:});
    
    plotit = p.Results.plotit;
    minspeed = p.Results.minspeed; 
    nBins = p.Results.nBins; 
    neurons = p.Results.neurons; 
    saveBool = p.Results.saveBool;
    shuffle = p.Results.shuffle; 
    noTreadmill = p.Results.noTreadmill;
    plotTrials = p.Results.plotTrials; 
    
    B = 1000;
    
    if isempty(neurons)
        neurons = getPlaceCells(md,.01);
    end
    if size(neurons,1) > size(neurons,2), neurons = neurons'; end 
    nNeurons = length(neurons);
    
%%
    currDir = pwd;
    cd(md.Location); 
    load('FinalOutput.mat','NumNeurons'); 
    nNeuronsTotal = NumNeurons; clear NumNeurons;
    
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
    load(fullfile(pwd,'Pos_align.mat'),'x_adj_cm','y_adj_cm','speed',...
        'time_interp','PSAbool','DFDTtrace');
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
    binocc = nan(size(X));
    parsed = postrials_treadmill(md);
    trials = parsed.summary(:,1)'; 
    nTrials = max(trials);
    
    keepgoing = true;
    bins = [1:0.001:nBins]';
    [F,raster] = deal(zeros(nTrials,nBins,nNeuronsTotal));
    curve = zeros(nNeuronsTotal,nBins);
    sigCurve = false(nNeuronsTotal,nBins); 
    smoothCurve = zeros(nNeuronsTotal,size(bins,1)); 
    
    [~,edges] = histcounts(X,nBins);
    
    for nn = neurons  
        if ~isnan(nn) && nn>0
            [raster(:,:,nn),F(:,:,nn),binocc] = ...
                MakeRasters(X,parsed,PSAbool(nn,:),DFDTtrace(nn,:),...
                isrunning,onTM,edges);

        curve(nn,:) = nanmean(raster(:,:,nn));
        smoothfit = fit([1:nBins]',curve(nn,:)','smoothingspline');
        smoothCurve(nn,:) = feval(smoothfit,bins)';  
        end
    end

%% Shuffle locations and re-compute tuning curves. 
    if shuffle
        shuffleCurve = nan(B,nBins,nNeuronsTotal);  %Shuffled curves. 
        shuffleCI = nan(2,nBins,nNeuronsTotal);     %Shuffled curve's confidence intervals.
        idxci = round([0.975;0.025].*B);            %Indices for confidence intervals.
        for nn = neurons
            for i=1:B
                %Make the raster. 
                shuffleRaster = permutePlace(raster(:,:,nn));
                
                %Take the mean across trials. 
                shuffleCurve(i,:,nn) = nanmean(shuffleRaster); 
            end
            %Sort then extract confidence intervals.
            temp = sort(shuffleCurve(:,:,nn)); 
            shuffleCI(:,:,nn) = temp(idxci,:); 
            
            %P-value for each data point. 
            p = sum(shuffleCurve(:,:,nn) > repmat(curve(nn,:),B,1))./B;
            
            sigCurve(nn,:) = p < 0.01; 
        end
    else 
        load(fullfile(md.Location),'SpatialTraces.mat','shuffleCurve',...
            'shuffleCI','sigCurve');
    end
          
    
% %% Get significance bins.
%     if calcSig
%         load('PlacefieldStats.mat','PFepochs','PFnHits'); 
%         [~,bestPF] = max(PFnHits,[],2);
% 
%         for thisNeuron = 1:nNeurons        
%             epochs = PFepochs{neurons(thisNeuron),bestPF(neurons(thisNeuron))};
%             nEpochs = size(epochs,1);
% 
%             for e = 1:nEpochs
%                 thisEpoch = epochs(e,1):epochs(e,2);
% 
%                 goodBins = binocc(thisEpoch);
%                 goodBins(isnan(goodBins)) = [];
% 
%                 sigCurve(neurons(thisNeuron),goodBins) = true;
%             end
% 
%         end
%     end

%% Plot
    if plotit
        thisNeuron = 1;
        while keepgoing
             
        shuffmean = mean(shuffleCurve(:,:,neurons(thisNeuron)));
        CImean = interp1(1:nBins,shuffmean,bins,'pchip'); 
        CIu = interp1(1:nBins,shuffleCI(1,:,neurons(thisNeuron)),bins,'pchip');
        CIl = interp1(1:nBins,shuffleCI(2,:,neurons(thisNeuron)),bins,'pchip');
        
        %Get spatial bins of significance.
        [SIGX,SIGY] = significance(bins,sigCurve(neurons(thisNeuron),:),...
                smoothCurve(neurons(thisNeuron),:),bins);
           
        f = figure(49); 
        f.Position = [550 180 360 565];
        subplot(2,1,1); 
            imagesc(raster(:,:,neurons(thisNeuron))); 
            colormap hot; caxis([0 1.2]);
            set(gca,'xtick',[],'ytick',[1,max(trials)],'linewidth',4,...
                'fontsize',12);
            
            if noTreadmill, xlim([22 80]); end
            
            title(['Neuron #',num2str(neurons(thisNeuron))],'fontsize',15);
    
            ylabel('Laps','fontsize',15); 
        subplot(2,1,2); hold on;                  
            if plotTrials, yyaxis left; end
            %Plot tuning curve. 
            plot(bins,smoothCurve(neurons(thisNeuron),:),...
                'color',[.58 .44 .86],'linewidth',5);
            
            %Plot confidence intervals. 
            plot(bins,CImean,'-b','linewidth',2);
            plot(bins,CIu,'--b',bins,CIl,'--b');
            
            %Plot significance markers. 
            Ylim = get(gca,'ylim');
            plot(SIGX,SIGY+Ylim(2)*0.1,'ro','linewidth',4);
            
            xlabel('Track position (cm)','fontsize',15);
            ylabel('Ca transient rate','fontsize',15);
            yLims = get(gca,'ylim');
            ylim([0 yLims(2)]);
            set(gca,'ycolor',[.58 .44 .86],'xtick',[1 80],'xticklabel',[0 200]);

            if noTreadmill
                xlim([22 80]);
                set(gca,'xtick',[22 80],'xticklabel',[0 140]);
                set(gca,'ycolor','k');
            end
            
            if plotTrials
                yyaxis right;
                plot(1:nBins,F(:,:,neurons(thisNeuron)),'-',...
                    'color',[.7 .7 .7 .2],'linewidth',2);
                set(gca,'ycolor',[.7 .7 .7]);
            end
            set(gca,'linewidth',4,'tickdir','out','fontsize',12);
            
            
        [keepgoing,thisNeuron] = scroll(thisNeuron,nNeurons,f);
        close all;
        
        end
    end
    
    cd(md.Location);
    if saveBool
        save('SpatialTraces.mat','raster','curve','smoothCurve',...
            'shuffleCurve','shuffleCI','sigCurve','nBins');
    end
    
    cd(currDir);
    
end

function [raster,F,binocc] = MakeRasters(X,parsed,PSAbool,trace,isrunning,onTM,edges)
    trials = parsed.summary(:,1)';
    nTrials = length(trials); 
    nBins = length(edges)-1;
    
    [F,raster] = deal(zeros(nTrials,nBins));
    for thisTrial = trials
        spkpos = X(parsed.trial==thisTrial & PSAbool & isrunning & ~onTM);

        [occ,~,binocc(parsed.trial==thisTrial)] = ...
            histcounts(X(parsed.trial==thisTrial),edges);
        binned = histcounts(spkpos,edges);

        binOccThisTrial = binocc(parsed.trial==thisTrial);
        traceThisTrial = trace(parsed.trial==thisTrial); 
        F(thisTrial,:) = accumarray(binOccThisTrial',...
            traceThisTrial',[nBins,1],@mean);
        raster(thisTrial,:) = binned./occ;
    end
end

function shuffled = permutePlace(X)
%shuffled = permuteTime(raster)
%
%   Performs a randomized card shuffle permutation for each row in raster. 
%
%   INPUT
%       raster: Trials x Time matrix of anything.
%
%   OUTPUT
%       shuffled: Same matrix, after shuffling elements in each row by a
%       random value. 
%
    
%% Shuffle. 
    [nLaps,nBins] = size(X);
    r = rem(randi([1,80],nLaps,1),nBins);
    wrapper = [X,X]; 
    shuffled = wrapper(bsxfun(@plus,bsxfun(@plus,nBins-r,0:nBins-1)*nLaps,(1:nLaps)'));
end

function [SIGX,SIGY] = significance(t,sigcurve,smoothedcurve,bins)
    %Get indices of significance.
    sigP = find(sigcurve);
    
    %Find the corresponding index in the smoothed bins.
    [~,inds] = ismember(sigP,bins);
    
    %Values to plot.
    SIGX = t(inds);
    SIGY = smoothedcurve(inds);
end