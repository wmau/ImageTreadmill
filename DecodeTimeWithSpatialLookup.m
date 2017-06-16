function DecodeTimeWithSpatialLookup(md,varargin)
%
%
%

%% Parse inputs. 
    path = md.Location;
    cd(path); 
    
    p = inputParser;
    p.addRequired('md'); 
    load(fullfile(path,'FinalOutput.mat'),'NumNeurons');
    p.addParameter('neurons',1:NumNeurons,@(x) isnumeric(x)); 
    p.addParameter('subsample',true,@(x) islogical(x)); 
    
    p.parse(md,varargin{:});
    
    neurons = p.Results.neurons;
    subsample = p.Results.subsample; 
    B = 1000;
    purple = [.58 .44 .86];
    teal = [0 .5 .5];

%% Bin position on treadmill.
    load(fullfile(path,'Pos_align.mat'),'x_adj_cm','y_adj_cm');
    x = x_adj_cm; y = y_adj_cm; 
    clear x_adj_cm y_adj_cm;
    load(fullfile(path,'TimeCells.mat'),'TodayTreadmillLog','T');
    
    %Get the treadmill indices and throw them into a vector. 
    good = [];
    [inds,nRuns] = TrimTrdmllInds(TodayTreadmillLog,T);
    for r=1:nRuns
        good = [good, inds(r,1):inds(r,2)];
    end
    
    %Set up variables for binning. 
    lims = [min(x(good)), max(x(good)); min(y(good)), max(y(good))];    %Set limits to treadmill boundaries.
    isrunning = true(size(x));                                          %Effectively always running.
    cmperbin = 1;
        
    %Bin position on the treadmill.
    [OccMap,~,~,~,xBin,yBin] = MakeOccMap(x,y,lims,isrunning,isrunning,cmperbin);

    %Set the limits in X and Y so that we can get linear indices.
    treadmillDims = size(OccMap);
    bin = zeros(size(x));               %Most of these entries are empty, but we need them to match treadmill indices. 
    bin(good) = sub2ind(treadmillDims,xBin(good),yBin(good));
    
%% Make lookup table.
    load(fullfile(path,'TreadmillTraces.mat'),'RawTrdmll'); 
    [~,nTimeBins,~] = size(RawTrdmll);
    nNeurons = length(neurons); 
    t = linspace(0,T,nTimeBins);
    
    %Get position in terms of spatial bin. 
    position = nan(nRuns,nTimeBins);        %Matrix containing bin per lap per frame.
    for r=1:nRuns
        position(r,:) = bin(inds(r,1):inds(r,2));
    end
    nBins = max(position(:));
    
    %Compute lookup table by taking the mean of fluorescence at that bin.
    lookup = nan(NumNeurons,nBins);           %NxT matrix containing lookup values. 
    for n=1:nNeurons
        traces = RawTrdmll(:,:,neurons(n));
        lookup(neurons(n),:) = accumarray(position(:),traces(:),[nBins,1],@mean)';
    end
    
%% Decode
    keepgoing = true; 
    n = 1;
    
    %Get the approximate fluorescence value at that spatial location.
    if subsample 
        pLapSample = 0.5;
        nLapSample = round(pLapSample*nRuns);
        
        [decoded,observed] = deal(nan(B,nTimeBins,NumNeurons));
        while keepgoing
            f = figure('Position',[680 560 265 420]); hold on;
            for i=1:B
                laps = randsample(1:nRuns,nLapSample);
                
                
                predictedTraces = reshape(lookup(neurons(n),position(laps,:)),...
                    [nLapSample,nTimeBins]);
                decoded(i,:,neurons(n)) = mean(predictedTraces);
                realTraces = RawTrdmll(laps,:,neurons(n));
                observed(i,:,neurons(n)) = mean(realTraces);
            end
            
            %Sort bootstrap samples for CI then take the mean. 
            sortedDecode = sort(decoded(:,:,neurons(n)));
            sortedObserved = sort(observed(:,:,neurons(n)));
            decodeM = mean(decoded(:,:,neurons(n)))';
            observedM = mean(observed(:,:,neurons(n)))';
            
            %Confidence intervals. 
            decodeL = sortedDecode(round(0.025*B),:)'; 
            decodeL = decodeM - decodeL; 
            decodeU = sortedDecode(round(0.975*B),:)'; 
            decodeU = decodeU - decodeM; 
            
            %Confidence intervals for observed traces. 
            observedL = sortedObserved(round(0.025*B),:)';
            observedL = observedM - observedL; 
            observedU = sortedObserved(round(0.975*B),:)';
            observedU = observedU - observedM; 
            
            %Plot.
            h = boundedline(t,decodeM,[decodeL,decodeU],...
                t,observedM,[observedL,observedU],'alpha'); 
            h(1).Color = purple;
            h(1).LineWidth = 2;
            h(2).Color = teal;
            h(2).LineWidth = 2; 
            set(gca,'tickdir','out','fontsize',12,'linewidth',4);
            xlabel('Time (s)','fontsize',15);
            ylabel('Fluorescence (A.U.)','fontsize',15);
            title(['Neuron #',num2str(neurons(n))],'fontsize',15);
            
            [keepgoing,n] = scroll(n,nNeurons,f); 
            close all;
            
        end
    else
        decoded = nan(nRuns,nTimeBins,nNeurons);
        for n=1:nNeurons
            for r=1:nRuns
                decoded(r,:,n) = lookup(neurons(n),position(r,:));
            end
        end
    end
end