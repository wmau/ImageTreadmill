function map = msMatchMultiSessionCells(mds,neurons)
%map = msMatchMultiSessionCells(mds,neurons)
%   
%   Gets a registration map of all neurons specified in the cell array
%   "neurons". 

%% 
    nSessions = length(mds);
    
    %Get the batch session map.
    mapMD = getMapMD(mds);
    cd(mapMD.Location);
    load('batch_session_map.mat');
    MAP = batch_session_map.map(:,2:end);
    
    %Preallocate.
    mapCols = zeros(nSessions,1);
    rows = [];    
    for s=1:nSessions
        cd(mds(s).Location);
    
        [~,mapRows,mapCols(s)] = msMatchCells(mds(s),neurons{s},false);
        rows = [rows; mapRows];
    end
    
    %No copies. 
    rows = unique(rows);
    
    %Get only the neurons and sessions specified. 
    map = MAP(rows,mapCols);
end