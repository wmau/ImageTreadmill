function [DS,sig] = DecodePlaceWithTemporalLookup(md,varargin)
%
%
%

%%
    path = md.Location;
    cd(path); 
    
    p = inputParser; 
    p.addRequired('md');
    load(fullfile(path,'FinalOutput.mat'),'NumNeurons');
    p.addParameter('neurons',1:NumNeurons,@(x) isnumeric(x)); 
    p.addParameter('neuralActivity','transient',@(x) ischar(x)); 
    p.addParameter('plotit',true,@(x) islogical(x)); 
    p.addParameter('noTreadmill',true,@(x) islogical(x));
    
    p.parse(md,varargin{:});
    
    neurons = p.Results.neurons;
    if size(neurons,1) > size(neurons,2), neurons = neurons'; end
    neuralActivity = p.Results.neuralActivity; 
    plotit = p.Results.plotit;
    noTreadmill = p.Results.noTreadmill;
    
    B = 100;
    purple = [.58 .44 .86];
    teal = [0 .5 .5];
    nNeurons = length(neurons);
    
%% 
    %Linearize position and bin it. Also extract trial timestamps.
    [~,~,~,~,X,parsed] = LinearizedPF_raster(md,'plotit',false,'neurons',[]);
    
    %Determine the number of trials and longest trial. 
    trialDurs = tabulate(parsed.trial); 
    nTrials = max(parsed.trial);
    longestTrialDuration = max(trialDurs(:,2));
    nSpatialBins = max(X); 
%%
    %Load fluorescence data. 
    switch neuralActivity
        case 'trace'
            load(fullfile(path,'Pos_align.mat'),'RawTrace'); 
            raster = RawTrace;
        case 'transient'
            load(fullfile(path,'Pos_align.mat'),'PSAbool'); 
            raster = PSAbool;
    end
    
    %Preallocate.
    position = nan(nTrials,longestTrialDuration); 
    trace = nan(nTrials,longestTrialDuration,NumNeurons);
    lookup = nan(NumNeurons,longestTrialDuration);
    
    %Get position and trace information. 
    for r=1:nTrials
        thisTrial = find(parsed.trial==r);
        dur = length(thisTrial);
        
        %Get position and throw it into the matrix. 
        position(r,1:dur) = X(thisTrial);
        
        for n=1:nNeurons
            trace(r,1:dur,neurons(n)) = raster(neurons(n),thisTrial);
        end
    end
    
%% Make lookup table. 
    %Simply average across trials. 
    for n=1:nNeurons
        lookup(neurons(n),:) = nanmean(trace(:,:,neurons(n)));
    end
    
%% 
    %Preallocate and some parameters. 
    pLapSample = 0.5;                       %Proportion of sampled laps during bootstrapping.
    nLapSample = round(pLapSample*nTrials); %Number of laps sampled during bootstrapping.
    
    if noTreadmill
        x = 22:nSpatialBins;
    else
        x = 1:nSpatialBins;                     %Spatial bins for plotting. 
    end
    
    if noTreadmill
        [decodedCurve,observedCurve] = deal(nan(NumNeurons,nSpatialBins-22+1));
        [dIterations,oIterations] = deal(nan(B,nSpatialBins-22+1,NumNeurons));
    else
        [decodedCurve,observedCurve] = deal(nan(NumNeurons,nSpatialBins));  
        [dIterations,oIterations] = deal(nan(B,nSpatialBins,NumNeurons));
    end
    sig = nan(NumNeurons,1); 
    keepgoing = true;   
    n = 1;
    
    if plotit && nNeurons > 0
        while keepgoing
            nn = neurons(n);
            f = figure('Position',[680 560 265 420]); hold on;
            
            %Make the curve. 
            [dC,oC,dI,oI,sig(nn)] = ...
                MakeCurve(position,lookup,trace,nn,nLapSample,x,B);
            
            if noTreadmill
                decodedCurve(nn,:) = dC(22:end);
                observedCurve(nn,:) = oC(22:end);
                dIterations(:,:,nn) = dI(:,22:end); 
                oIterations(:,:,nn) = oI(:,22:end); 
            else 
                decodedCurve(nn,:) = dC;
                observedCurve(nn,:) = oC;
                dIterations(nn,:) = dI;
                oIterations(nn,:) = oI; 
            end

            
            %Sort the matrices for making CIs. 
            sortedDecode = sort(dIterations(:,:,nn)); 
            sortedObserved = sort(oIterations(:,:,nn)); 

            %Get confidence intervals. 
            decodeCI = [    decodedCurve(nn,:)' - sortedDecode(round(0.025*B),:)',...
                            sortedDecode(round(0.975*B),:)' - decodedCurve(nn,:)'];
            observedCI = [  observedCurve(nn,:)' - sortedObserved(round(0.025*B),:)',...
                            sortedObserved(round(0.975*B),:)' - observedCurve(nn,:)']; 

            %Plot. 
            [h,p] = boundedline(x,decodedCurve(nn,:),decodeCI,...
                x,observedCurve(nn,:),observedCI,...
                'nan','gap','alpha'); 
            h(1).Color = teal; 
            p(1).FaceColor = teal;
            h(2).Color = purple; 
            p(2).FaceColor = purple;
            [h.LineWidth] = deal(2); 
            if noTreadmill
                set(gca,'xtick',[22 80],'xticklabel',[0 140]);
            else
                set(gca,'xtick',[1 80],'xticklabel',[0 200]); 
            end
            set(gca,'tickdir','out','fontsize',12,'linewidth',4);
            xlabel('Position (cm)','fontsize',15); 
            switch neuralActivity
                case 'trace'
                    ylabel('Fluorescence (A.U.)','fontsize',15);
                case 'transient'
                    ylabel('Ca^{2+} transient rate','fontsize',15);
            end
            title(['Neuron #',num2str(nn)],'fontsize',15); 

            [keepgoing,n] = scroll(n,nNeurons,f); 
            close all;
        end
    elseif nNeurons > 0
        for nn=neurons
            [dC,oC,~,~,sig(nn)] = MakeCurve(position,lookup,trace,nn,nLapSample,x,B);
            
            if noTreadmill
                decodedCurve(nn,:) = dC(22:end);
                observedCurve(nn,:) = oC(22:end);
            else 
                decodedCurve(nn,:) = dC;
                observedCurve(nn,:) = oC;
            end
        end
    end

%% Normalize areas under the curve then find difference. 
    [dNorm,oNorm] = deal(nan(size(decodedCurve)));
    DS = nan(NumNeurons,1);
    if nNeurons > 0
        for nn=neurons
            %Non-negative values only.
            if any(decodedCurve(nn,:)<0)
                dNorm(nn,:) = decodedCurve(nn,:) + abs(min(decodedCurve(nn,:)));
            else 
                dNorm(nn,:) = decodedCurve(nn,:);
            end

            if any(observedCurve(nn,:)<0)
                oNorm(nn,:) = observedCurve(nn,:) + abs(min(observedCurve(nn,:)));
            else
                oNorm(nn,:) = observedCurve(nn,:);
            end

            %Normalize so area under the curve = 1.
            dNorm(nn,:) = dNorm(nn,:)./trapz(x,dNorm(nn,:));
            oNorm(nn,:) = oNorm(nn,:)./trapz(x,oNorm(nn,:)); 

            %Find difference. 
            DS(nn) = trapz(x,abs(dNorm(nn,:)-oNorm(nn,:)));     
        end
    end

end

function [decodedCurve,observedCurve,dIterations,oIterations,sig] = ....
    MakeCurve(position,lookup,trace,neuron,nLaps,x,B)

    nSpatialBins = max(position(:));
    nTrials = size(position,1);
    [dIterations,oIterations] = deal(nan(B,nSpatialBins));
   
    for i=1:B
        %Randomly sample some laps. 
        laps = randsample(1:nTrials,nLaps);
        subset = position(laps,:);

        %For each bin, get timestamps the mouse was in that position,
        %then consult the lookup table to find average fluorescence
        %response at that time.
        [onMaze,binMatch] = ismember(subset,x); 
        [~,t] = find(onMaze); 
        binMatch = binMatch(binMatch~=0); 

        curve = accumarray(binMatch,lookup(neuron,t),[nSpatialBins,1],@nanmean);
        dIterations(i,:) = curve';

        %Only consider non-nans.
        tempPos = position(laps,:);
        tempTrace = trace(laps,:,neuron); 
        good = ~isnan(tempPos); 

        %Make empirical tuning curve. 
        oIterations(i,:) = accumarray(tempPos(good),tempTrace(good),[nSpatialBins 1],@nanmean); 
    end
    
    decodedCurve = nanmean(dIterations); 
    observedCurve = nanmean(oIterations); 
    
    temp = repmat(observedCurve,[B,1]); 
    p = sum(temp<dIterations)./B; 
    sig = any(p<0.05);
end