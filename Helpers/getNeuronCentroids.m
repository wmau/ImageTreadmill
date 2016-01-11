function centroids = getNeuronCentroids(animal,date,session)
%centroids = getNeuronCentroids(animal,date,session)
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
    ChangeDirectory(animal,date,session);
    
    %Load neuron masks. 
    try
        load(fullfile(pwd,'ProcOut.mat'),'NeuronImage');
    catch
        disp('ProcOut.mat not found. Run TENASPIS!'); 
    end
    
    %Get the centroid of the mask. 
    props = cellfun(@regionprops,NeuronImage); 
    temp = extractfield(props,'Centroid'); 
    
    %Every other element is an X or Y coordinate. 
    centroids = [temp(1:2:end)', temp(2:2:end)'];
   
end