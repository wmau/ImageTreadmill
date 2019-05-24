clear; 
clc;
loadMD;

sessions = [MD(292:303) MD(305:308)];   
n_sessions = length(sessions);

[Traces, Whole_Traces, Whole_Transients, Treadmill_Indices, Transients, ...
    Time_Cells, Significant, STBS] = deal(cell(4));
for session_number = 1:n_sessions
    % Change directory. 
    this_session = sessions(session_number); 
    cd(this_session.Location); 
    
    % Load relevant variables. 
    load(fullfile(this_session.Location, 'TimeCells.mat'), ...
        'TodayTreadmillLog', 'curves', ...
        'ratebylap', 'T');
    load(fullfile(this_session.Location, 'TreadmillTraces.mat'), ...
        'RawTrdmll');
    load(fullfile(this_session.Location, 'Pos_align.mat'), ...
        'PSAbool','RawTrace'); 
    n_neurons = size(PSAbool,1);   
    
    % Get time cells. 
    TimeCells = AcquireTimePlaceCells(this_session, 'timecells')'; 
    
    % Make raster. 
    inds = TrimTrdmllInds(TodayTreadmillLog, T);
    raster = nan(size(RawTrdmll)); 
    
    for this_neuron = 1:n_neurons
        raster(:,:,this_neuron) = buildRasterTrace(inds, ...
            PSAbool, this_neuron); 
    end
    
    % Get session trial bias scores.
    [skewness, ~, ~, sigCurve] = getAllSkewnesses(this_session);
    
    % Distribute. 
    Traces{session_number} = RawTrdmll; 
    Whole_Traces{session_number} = RawTrace;
    Whole_Transients{session_number} = PSAbool;
    Treadmill_Indices{session_number} = inds; 
    Transients{session_number} = raster;
    Time_Cells{session_number} = TimeCells; 
    Significant{session_number} = sigCurve; 
    STBS{session_number} = skewness';
end

% Get animal's registration cell map. 
Maps = cell(4,1);
c = 1;
for session_number = 1:4:n_sessions
    % Get sessions and load map.
    animal_sessions = sessions(session_number:session_number+3);
    dates = {animal_sessions.Date};
    session_numbers = [animal_sessions.Session];
    n_sessions = length(animal_sessions);
    
    map_MD = getMapMD(animal_sessions); 
    cd(map_MD.Location); 
    load('batch_session_map.mat'); 
    
    % Order map to be chronological. 
    regDates = {batch_session_map.session.Date};
    regSessions = [batch_session_map.session.Session];
    map =  batch_session_map.map(:,2:end); 
    
    mapCols = zeros(n_sessions,1);
    for i=1:n_sessions
        mapCols(i) = find(ismember(regDates,dates{i}) ...
            & ismember(regSessions,session_numbers(i)));
    end
    
    map = map(:, mapCols);
    Maps{c} = map; 
    
    c = c+1;
end

Traces = Traces';
Whole_Traces = Whole_Traces';
Whole_Transients = Whole_Transients';
Treadmill_Indices = Treadmill_Indices';
Transients = Transients';
Time_Cells = Time_Cells';
Significant = Significant';
STBS = STBS'; 
Maps = Maps';

save(fullfile('C:\Users\William Mau\Desktop\Data for Yue', ...
    'YueData.mat'), 'Traces', 'Whole_Traces', 'Whole_Transients',...
    'Treadmill_Indices', 'Transients', 'Time_Cells', 'Significant', ...
    'STBS', 'Maps');