function batch_tenaspis(MD)
%batch_tenaspis(MD)
%
%   Runs TENASPIS on multiple sessions. 
%
%   INPUT
%       MD: Master Directory entries. 
%

%% Change directory and run TENASPIS. 
    nSessions = length(MD); 
    
    for i=1:nSessions
        cd(fullfile(MD(i).Location,'ICmovie_smoothed-Objects'));
        infile = dir('*.h5');
        
        Tenaspis(infile,'animal_id',MD(i).Animal,'sess_date',MD(i).Date,'sess_num',MD(i).Session);
    end
    
end