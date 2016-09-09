function md = findMDfromGraphData(graphData)
%md = findMDfromGraphData(graphData)
%
%   graphData has the animal, date, and session embedded in its fields.
%   This function looks for the correct entry from MD that matches those
%   entries. 
%
%   INPUT
%       graphData: from MakeGraphv4.
%
%   OUTPUT
%       md: session entry. 

%% Find MD. 
    %Load MD and compile all the data. 
    loadMD;
    animals = {MD.Animal};
    dates = {MD.Date};
    sessions = [MD.Session];
    
    %Find entry that matches all criteria. 
    entry = find(strcmp(animals,graphData.Animal) & ...
        strcmp(dates,graphData.Date) & ...
        sessions == graphData.Session);
    
    %Output. 
    md = MD(entry); 
end