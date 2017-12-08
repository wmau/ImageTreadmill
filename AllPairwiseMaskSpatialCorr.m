function [d,r,matchMat] = AllPairwiseMaskSpatialCorr(base,reg,varargin)
%[d,r,matchMat] = AllPairwiseMaskSpatialCorr(base,reg)
%
%   To provide convincing evidence that we correctly segment and register
%   the same cells across days, we perform spatial correlations of all the
%   ROIs on one day to all the ROIs on another day. If the collection of
%   cell bodies is so dense that we would not be able to distinguish
%   between neighboring cells, then we would see a continuum of correlation
%   values that gradually decreases as a function of distance. However, if
%   we observe a cluster of values at high correlation, low distance that
%   is clearly distinct from the rest of the values, then we can assert
%   that the closest neurons across days are also most likely the same
%   neurons since they share the same shape.
%
%   INPUTS
%       base & reg: Base and registered sessions. 
%
%   OUTPUTS
%       d: Matrix of pairwise distances between centroids (in microns).
%       Rows = cells in base session. Columns = cells in registered
%       session.
%       
%       r: Matrix of pairwise ROI spatial corerlations (r). 
%
%       matchMat: Two-column vector specifying cells in base session
%       (column 1) and their index in the registered session (column 2).
%

%% Parse inputs.
    p = inputParser; 
    p.addRequired('base',@(x) isstruct(x)); 
    p.addRequired('reg',@(x) isstruct(x)); 
    p.addParameter('corr_vs_dist',true,@(x) islogical(x)); 
    p.addParameter('corr_vs_nextNearest',true,@(x) islogical(x)); 
    
    p.parse(base,reg,varargin{:});
    corr_vs_dist = p.Results.corr_vs_dist; 
    corr_vs_nextNearest = p.Results.corr_vs_nextNearest; 

%% Align the fields of view from the two sessions. 
    %Get registration information.
    regInfo = image_registerX(base.Animal,base.Date,base.Session,reg.Date,reg.Session);
    
    %Load base session's ROIs.
    cd(base.Location);
    load('FinalOutput.mat','NeuronImage','NeuronAvg'); 
    baseNeuronImage = NeuronImage; 
    
    %Make the average value of the ROI during calcium transient activity.
    baseNeuronAvg = MakeAvgROI(NeuronImage,NeuronAvg); 
    
    %Load registered session's ROIs. 
    cd(reg.Location);
    load('FinalOutput.mat','NeuronImage','NeuronAvg'); 
    
    %Registered session's ROIs after transformation. 
    regNeuronImage = cellfun(@(x) imwarp_quick(x,regInfo),NeuronImage,...
        'UniformOutput',false); 
    
    %Make the average value of the ROI during calcium transient activity.
    regNeuronAvg = MakeAvgROI(NeuronImage,NeuronAvg);
    
    %Transform to align to base session.
    regNeuronAvg = cellfun(@(x) imwarp_quick(x,regInfo),regNeuronAvg,...
        'UniformOutput',false); 

%% Get cell centroids.
    %Preallocate. The first and second cell arrays in centroids correspond
    %to the base and registered sessions respectively. 
    centroids = cell(1,2);
    
    %Get centroids of the ROIs. These are necessary to center them. 
    props = cellfun(@(x) regionprops(x,'Centroid'),baseNeuronImage);
    temp = extractfield(props,'Centroid'); 
    
    %Centroids for base session.
    centroids{1} = [temp(1:2:end)', temp(2:2:end)'];
    
    %Get centroids of ROIs. 
    props = cellfun(@(x) regionprops(x,'Centroid'),regNeuronImage); 
    temp = extractfield(props,'Centroid'); 
    
    %Centroids for registered session.
    centroids{2} = [temp(1:2:end)', temp(2:2:end)'];
    
%% Make window around neuron ROI to correlate with other ROIs. 
    %Preallocate. 
    nNeuronsBase = length(baseNeuronImage);     %Number of neurons in base session,
    nNeuronsReg = length(regNeuronImage);       %Number of neurons in registered session. 
    [xWindow,yWindow,windowedBase] = deal(cell(nNeuronsBase,1));
    dims = size(baseNeuronImage{1}); 
    
    %Make a 16-pixel window around the neuron to correlate with other masks
    for n=1:nNeuronsBase
        
        %Columns in the image.
        xLowerBound = round(centroids{1}(n,1)-8);   %Leftmost boundary.
        xUpperBound = xLowerBound + 16;             %Rightmost boundary. 
        
        %To handle cases on the edge. 
        if xLowerBound < 1, xLowerBound = 1; end                %Ran out of space on the left?
        if xUpperBound > dims(2), xUpperBound = dims(2); end    %Ran out of space on the right?
        
        %Make the window. 
        xWindow{n} = xLowerBound:xUpperBound;
        
        %Rows in the image. 
        yLowerBound = round(centroids{1}(n,2)-8);   %Bottommost boundary.
        yUpperBound = yLowerBound + 16;             %Topmost boundary. 
        
        %To handle cases on the edge. 
        if yLowerBound < 1, yLowerBound = 1; end                %Ran out of space at the bottom? 
        if yUpperBound > dims(1), yUpperBound = dims(1); end    %Ran out of space at the top?
        
        %Make the window. 
        yWindow{n} = yLowerBound:yUpperBound; 
        
        %Make the window around the neuron. 
        windowedBase{n} = baseNeuronAvg{n}(yWindow{n},xWindow{n});
    end

%% Get distances between ROI centroids. 
    %Preallocate. 
    [r,d] = deal(zeros(nNeuronsBase,nNeuronsReg));
    
    %Get distances. 
    for n1=1:nNeuronsBase
        x1 = centroids{1}(n1,1);        %x-coordinate of base session neuron centroid.
        y1 = centroids{1}(n1,2);        %y-coordinate of base session neuron centroid. 
              
        for n2=1:nNeuronsReg
            x2 = centroids{2}(n2,1);    %x-coordinate of registered session neuron centroid.
            y2 = centroids{2}(n2,2);    %y-coordinate of registered session neuron centroid. 
        
            %Get distance between centroids. 
            conversion = 1.1;           %1.1 microns per pixel. 
            d(n1,n2) = conversion*sqrt((x2-x1)^2 + (y2-y1)^2);
        end
    end

%% Correlate neuron masks. 
    for n1=1:nNeuronsBase  
        %Get window used for the base neuron ROI. 
        tempY = yWindow{n1}; 
        tempX = xWindow{n1};
        
        for n2=1:nNeuronsReg
            %Use the same window as the base neuron's. 
            windowedReg = regNeuronAvg{n2}(tempY,tempX);
            
            %Spatial correlation. 
            r(n1,n2) = corr(windowedBase{n1}(:),windowedReg(:));
        end
    end
    
%% Get real neuron mappings.
    matchMat = msMatchCells([base,reg],1:nNeuronsBase,true);
    
%% Plot correlations vs. distances. 
    %Flatten. 
    distances = d(:);
    correlations = r(:);
 
    %Get the distances and correlation values from the cells that we
    %matched across days. 
    nMatches = size(matchMat,1);
    [matchedDs,matchedRs,nextNearestRs] = deal(zeros(nMatches,1));
    for n=1:nMatches
        matchedDs(n) = d(matchMat(n,1),matchMat(n,2)); 
        matchedRs(n) = r(matchMat(n,1),matchMat(n,2)); 
    end
           
    %Plotting ROI spatial correlations against centroid distance. 
    if corr_vs_dist
        %Only plot cell pairs across the two sessions that are under 20 microns
        %away.
        thresh = 20; 
        nearby = d<thresh;

        %Plot distance vs. correlation for all cell pairs across days (whose
        %distance is under threshold).
        figure;
        scatter(distances(nearby),correlations(nearby),'.'); 
        hold on;

        %Plot the matched cell pairs in green. 
        scatter(matchedDs,matchedRs,'g.');
        xlabel('Centroid distance (microns)');
        ylabel('ROI spatial correlation (r)');
        axis square;
        set(gca,'tickdir','out');
    end
   
%% Plot correlations of next nearest neighbor. 
    if corr_vs_nextNearest
        for n=1:nMatches
            sortedDistances = sort(d(matchMat(n,1),:)); 
            nextNearest = find(d(matchMat(n,1),:) == sortedDistances(2));
            
            nextNearestRs(n) = r(matchMat(n,1),nextNearest);
        end
        
        figure;
        scatter(matchedRs,nextNearestRs,'.');
    end
end