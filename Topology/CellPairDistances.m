function d = CellPairDistances(md,n1,n2,varargin)
%d = CellPairDistances(md,n1,n2,varargin)
%
%   Calculates Euclidean distance between ROI centroids of neurons listed
%   in n1 and n2. 
%
%   INPUTS
%       md: session entry.
%
%       n1 and n2: two vectors of neurons. 
%
%       NAME,VALUE:
%           centroids: Nx2 matrix of centroids, from getNeuronCentroids(),
%           useful if you are running CellPairDistances() more than once
%           since getNeuronCentroids is computationally taxing. 
%

%% Parse inputs. 

    p = inputParser; 
    p.addRequired('md',@(x) isstruct(x)); 
    p.addRequired('n1',@(x) isnumeric(x));
    p.addRequired('n2',@(x) isnumeric(x)); 
    p.addParameter('centroids',[],@(x) isnumeric(x)); 

    p.parse(md,n1,n2,varargin{:});
    
    centroids = p.Results.centroids;
%% Get cell pair distances. 
    cd(md.Location);
    
    %Make it a column vector. 
    if size(n1,2) < size(n1,1), n1 = n1';end
    if size(n2,2) < size(n2,1), n2 = n2';end
    
    neuronList = [n1,n2]; 
    if isempty(centroids)
        centroids = getNeuronCentroids(md,'neurons',neuronList);
    end
    
    %Neurons in each vector. 
    nNeurons1 = length(n1); 
    nNeurons2 = length(n2); 
    
    %Preallocate. 
    d = nan(nNeurons1,nNeurons2); 
    for i=1:nNeurons1
        %Centroids for cell 1.
        x1 = centroids(n1(i),1); 
        y1 = centroids(n1(i),2); 
        
        for j=1:nNeurons2
            %Centroid for cell 2.
            x2 = centroids(n2(j),1); 
            y2 = centroids(n2(j),2); 
            
            %Anatomical distance. 
            d(i,j) = sqrt((x2-x1)^2 + (y2-y1)^2);
            
        end
    end
   
end