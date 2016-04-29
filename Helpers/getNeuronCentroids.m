function centroids = getNeuronCentroids(MD,varargin)
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
    animal = MD.Animal;
    date = MD.Date;
    session = MD.Session;
    ChangeDirectory(animal,date,session);
    
    %Minimum transient length input into TENASPIS. 
    if nargin>3
        min_trans_length = varargin{1}; 
    end

    %Load neuron masks. 
    try
        if exist('min_trans_length','var')
            load(fullfile(pwd,['ProcOut_minlength_',num2str(min_trans_length),'.mat']),'NeuronImage');
        else
            load(fullfile(pwd,'ProcOut.mat'),'NeuronImage');
        end
    catch
        disp('ProcOut.mat not found. Run TENASPIS!'); 
    end
    
    %Get the centroid of the mask. 
    props = cellfun(@regionprops,NeuronImage); 
    temp = extractfield(props,'Centroid'); 
    
    %Every other element is an X or Y coordinate. 
    centroids = [temp(1:2:end)', temp(2:2:end)'];
   
end