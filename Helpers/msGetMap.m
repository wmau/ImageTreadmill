function matchedCells = msGetMap(mds)
%matchedCells = msGetMap(mds)
%
%   Simply gets the map of all neurons in each session in mds. 
%

%% Main function.
    nSessions = length(mds); 
    
    %Make a vector listing each cell number. 
    neurons = cell(1,nSessions); 
    for s=1:nSessions
        cd(mds(s).Location);
        load('FinalOutput.mat','NumNeurons'); 
        
        neurons{s} = 1:NumNeurons; 
    end
    
    %Get the map. 
    matchedCells = msMatchMultiSessionCells(mds,neurons);
end