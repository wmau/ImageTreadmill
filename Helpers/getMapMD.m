function mapMD = getMapMD(md)
%mapMD = getMapMD(md)
%
%   Gets the batch_session_map for the sessions indicated. These should
%   already go together. 
%
%   INPUT
%       md: session entry. These should be from one animal and the
%       registered sessions should already have been done. 
%
%   OUTPUT
%       mapMD: session entry containing the batch_session_map. 
%

%% Main body. 
    animals =   {md.Animal}; 
    dates =     {md.Date};
    sessions =  [md.Session]; 
    
    loadMD;
    %Insert batch_session_map directory here for additional animals.
    mds =      [MD(215);        %1 - G45
                MD(253);        %2 - G48
                MD(291);        %3 - Polaris
                MD(275)];       %4 - Bellatrix
        
    %Get registration data. 
    nAnimals = length(mds); 
    regData = cell(1,nAnimals); 
    animalStrMatch = false(1,nAnimals);
    for a=1:nAnimals
        cd(mds(a).Location);
        load('batch_session_map.mat');
        regData{a} = batch_session_map.session;
        
        %Do all animals in md match? Which mds entry?
        animalStrMatch(a) = all(strcmp(animals,mds(a).Animal));
    end
   
    %Compile the registration data. 
    regAnimals =    {regData{animalStrMatch}.Animal};
    regDates =      {regData{animalStrMatch}.Date};
    regSessions =   [regData{animalStrMatch}.Session];
    
    %Make sure all the entries in md match those in a batch_session_map.
    good =  all(ismember(animals,regAnimals) & ...
            ismember(dates,regDates) & ...
            ismember(sessions,regSessions));

    %If no good, error. 
    assert(good,'Error - at least one entry does not match batch_session_map'); 
    
    %Otherwise, spit out the session entry containing the
    %batch_session_map.
    mapMD = mds(animalStrMatch);
end