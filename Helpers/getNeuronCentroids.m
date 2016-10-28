function centroids = getNeuronCentroids(MD,varargin)
%centroids = getNeuronCentroids(MD,varargin)
%
%   INPUTS
%       animal: Name of the mouse (e.g., GCamp6f_45_treadmill).
%
%       date: e.g., 11_20_2015
%
%       session: Session number.
%
%   OUTPUT
%       centroids: Nx2 matrix where N is the number of neurons. First
%       column is the X coordinate. Second column is the Y coordinate. 
%

%% Obtain neuron centroids.
    cd(MD.Location); 
 
    %Load neuron masks. 
    load(fullfile(pwd,'FinalOutput.mat'),'NeuronImage');
    nNeurons = length(NeuronImage);
    
    %Parse inputs.
    p = inputParser;
    p.addRequired('MD',@(x) isstruct(x));
    p.addParameter('neurons',1:nNeurons,@(x) isnumeric(x));
    
    p.parse(MD,varargin{:});
    neurons = p.Results.neurons;
    
    %Get the centroid of the mask. 
    props = cellfun(@(x) regionprops(x,'Centroid'),NeuronImage(neurons)); 
    temp = extractfield(props,'Centroid'); 
    
    %Every other element is an X or Y coordinate. 
    centroids = nan(nNeurons,2);
    centroids(neurons,:) = [temp(1:2:end)', temp(2:2:end)'];
   
end