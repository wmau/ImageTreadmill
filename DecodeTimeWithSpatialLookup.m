function [DS] = DecodeTimeWithSpatialLookup(md,varargin)
%[DS] = DecodeTimeWithSpatialLookup(md,varargin)
%
%   

%% Parse inputs. 
    path = md.Location;
    cd(path); 
    
    p = inputParser;
    p.addRequired('md'); 
    load(fullfile(path,'FinalOutput.mat'),'NumNeurons');
    p.addParameter('neurons',1:NumNeurons,@(x) isnumeric(x)); 
    p.addParameter('plotit',true,@(x) islogical(x));
    
    p.parse(md,varargin{:});
    
    neurons = p.Results.neurons;
    if size(neurons,1) > size(neurons,2), neurons = neurons'; end
    plotit = p.Results.plotit;
    B = 100;
    purple = [.58 .44 .86];
    teal = [0 .5 .5];
    nNeurons = length(neurons); 

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
    [dIterations,oIterations] = deal(nan(B,nTimeBins,NumNeurons)); 
    [decodedCurve,observedCurve] = deal(nan(NumNeurons,nTimeBins));
    pLapSample = 0.5;
    nLapSample = round(pLapSample*nRuns);
    
    %Get the approximate fluorescence value at that spatial location.
    if plotit && nNeurons > 0
        while keepgoing                
            nn = neurons(n);
            
            %Make decoded and observed curves. 
            [decodedCurve(nn,:),observedCurve(nn,:),dIterations(:,:,nn),oIterations(:,:,nn)] = ....
                MakeCurve(position,lookup,RawTrdmll,nn,nLapSample,B);
            
            %Sort bootstrap samples for CI then take the mean. 
            sortedDecode = sort(dIterations(:,:,nn));
            sortedObserved = sort(oIterations(:,:,nn));

            %Confidence intervals. 
            decodeCI = [decodedCurve(nn,:)' - sortedDecode(round(0.025*B),:)',...
                        sortedDecode(round(0.975*B),:)' - decodedCurve(nn,:)']; 

            %Confidence intervals for observed traces. 
            observedCI = [   observedCurve(nn,:)' - sortedObserved(round(0.025*B),:)',...
                            sortedObserved(round(0.975*B),:)' - observedCurve(nn,:)']; 

            %Plot.
            f = figure('Position',[680 560 265 420]); hold on;
            [h,p] = boundedline(t,decodedCurve(nn,:),decodeCI,...
                t,observedCurve(nn,:),observedCI,'alpha'); 
            h(1).Color = purple;
            p(1).FaceColor = purple;
            h(1).LineWidth = 2;
            h(2).Color = teal;
            p(2).FaceColor = teal;
            h(2).LineWidth = 2; 
            set(gca,'tickdir','out','fontsize',12,'linewidth',4);
            xlabel('Time (s)','fontsize',15);
            ylabel('Fluorescence (A.U.)','fontsize',15);
            title(['Neuron #',num2str(neurons(n))],'fontsize',15);

            [keepgoing,n] = scroll(n,nNeurons,f); 
            close all;

        end
    elseif nNeurons > 0
        for nn=neurons
            [decodedCurve(nn,:),observedCurve(nn,:),dIterations(:,:,nn),oIterations(:,:,nn)] = ....
                MakeCurve(position,lookup,RawTrdmll,nn,nLapSample,B);
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
            dNorm(nn,:) = dNorm(nn,:)./trapz(t,dNorm(nn,:));
            oNorm(nn,:) = oNorm(nn,:)./trapz(t,oNorm(nn,:)); 

            %Find difference. 
            DS(nn) = trapz(t,abs(dNorm(nn,:)-oNorm(nn,:)));    
        end
    end
    
end

function [decodedCurve,observedCurve,dIterations,oIterations] = ....
    MakeCurve(position,lookup,traces,neuron,nLapSample,B)
%% Make the decoded and observed curve based on a lookup table. 

    %Dimensions.
    [nRuns,nTimeBins] = size(position);
    
    %Preallocate. 
    [dIterations,oIterations] = deal(nan(B,nTimeBins));    
    
    %Randomly sample trials.
    for i=1:B
        laps = randsample(1:nRuns,nLapSample);

        %Extract fluorescence values based on position on the treadmil.
        try
        predictedTraces = reshape(lookup(neuron,position(laps,:)),...
            [nLapSample,nTimeBins]);
        catch, keyboard; end
        
        %Tuning curve by taking the mean of predicted traces per lap. 
        dIterations(i,:) = mean(predictedTraces);
        
        %Get real traces then take the mean to get actual tuning curve. 
        realTraces = traces(laps,:,neuron);
        oIterations(i,:) = mean(realTraces);
    end
    
    %Take the mean over all the iterations of tuning curves. 
    decodedCurve = mean(dIterations)';
    observedCurve = mean(oIterations)';
end