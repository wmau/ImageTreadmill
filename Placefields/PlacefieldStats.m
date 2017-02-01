function PlacefieldStats(md)
%PlacefieldStats(md)
%
%   Calculates basic properties of place fields such as their regional
%   area and the percentage of 'hits' it gets as a proportion of passes. 
%   
%   INPUT
%       md: session to analyze.
%
%   OUTPUTS
%       Each of these are matrices or cell arrays of size NxP (N=# neurons,
%       P=# place fields).
%
%       PFpcthits: percentage of epochs that contained at least one
%       activation of the place cell in that place field.
%
%       PFnHits: number of times that cell was active in that place field.
%
%       PFnEpochs: number of passes through that place field. Each
%       traversal was only counted once (i.e., we are not counting all the
%       frames that the mouse was in that place field).
%
%       PFepochs: cell array containing all the epochs (start = 1st column,
%       end = 2nd column) where mouse passed through that place field.
%
%       PFcentroids: centroids of place fields.
%
%       PFpixels: pixels (in linear index form) of place fields.
%
%       PFarea: area of place field.
%
%       bestPF: vector indexing the column corresponding to the place field
%       with the highest activation.
%
%% Set up.
    cd(md.Location);
    load('Placefields.mat','TMap_gauss','xBin','yBin','isrunning');
    try
        load('Pos_align.mat','PSAbool');
    catch
        load('FinalOutput.mat','PSAbool');
        [~,~,~,PSAbool] = AlignImagingToTracking(md.Pix2CM,PSAbool,0);
    end
    PSAbool = PSAbool(:,isrunning);
    
%% Get basic properties of the placefields
    nNeurons = length(TMap_gauss);
    cc = cell(1,nNeurons);
    PFprops = cell(1,nNeurons);
    for n=1:nNeurons
        %Find peak.
        peak = max(TMap_gauss{n}(:));
        
        %Create binary image where anything above half the peak is 1.
        binImage = TMap_gauss{n} > peak/2;
        
        %Get place field blobs and basic properties.
        cc{n} = bwconncomp(binImage);
        PFprops{n} = regionprops(cc{n},'area','centroid'); 
    end
    
    %Get the number of place fields each cell has. 
    nPFs = cellfun(@(x) x.NumObjects,cc);
    maxNPFs = max(nPFs);
    
%% Get epochs of place field traversal.
    %Convert to linear indices.
    linInd = sub2ind(size(TMap_gauss{1}),xBin,yBin);
    
    %Preallocate a lot of shit. 
    PFpixels = cell(nNeurons,maxNPFs);
    PFarea = nan(nNeurons,maxNPFs);
    PFcentroids = cell(nNeurons,maxNPFs);
    PFepochs = cell(nNeurons,maxNPFs);
    PFnEpochs = zeros(nNeurons,maxNPFs);
    PFactive = cell(nNeurons,maxNPFs);
    bestPF = ones(nNeurons,1);
    
    for n=1:nNeurons
        %Compile place field pixels, area, and centroids.
        PFpixels(n,1:nPFs(n)) = cc{n}.PixelIdxList(:)';
        PFarea(n,1:nPFs(n)) = [PFprops{n}.Area];
        PFcentroids(n,1:nPFs(n)) = {PFprops{n}.Centroid};
               
        for p=1:nPFs(n)
            %For each place field, find when the mouse was in it.
            inPF = ismember(linInd,PFpixels{n,p});
            
            %Get traversal indices.
            PFepochs{n,p} = NP_FindSupraThresholdEpochs(inPF,eps,0);
            PFnEpochs(n,p) = size(PFepochs{n,p},1);
            
            PFactive{n,p} = zeros(PFnEpochs(n,p),1);
            for epoch=1:PFnEpochs(n,p)
                %Start and stop indices for traversal.
                s = PFepochs{n,p}(epoch,1);
                e = PFepochs{n,p}(epoch,2);
                
                %Get activations during traversal epochs.
                PFactive{n,p}(epoch) = any(PSAbool(n,s:e));                             
            end
        end
        
        %Get peak activation.
        [~,peakPix] = max(TMap_gauss{n}(:));
        
        %If there's a place field...
        if ~all(cellfun('isempty',PFpixels(n,:)))
            %Find each place field it corresponds to.
            bestPF(n) = find(cellfun(@(x) ismember(peakPix,x),PFpixels(n,:)));
        end
    end
    
    %Number and percentage hits.
    PFnHits = cellfun(@sum,PFactive);
    PFpcthits = PFnHits./PFnEpochs;
     
    save('PlacefieldStats.mat','PFpcthits','PFnHits','PFnEpochs','PFepochs',...
        'PFcentroids','PFpixels','PFarea','bestPF','-v7.3');
end