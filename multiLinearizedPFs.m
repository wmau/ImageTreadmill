function multiLinearizedPFs(mapMD,base,comp)
%multiLinearizedPFs(mapMD,base,comp)
%
%   Plots the place fields of place cells in 'base' session then the
%   spatial responses of those same cells over the sessions specified in
%   'comp'. Order determined by the peak response of 'base'. 
%
%   INPUTS
%       mapMD: MD entry where batch_session_map lives. 
%
%       base: The function will look for place cells in this entry and
%       organize all plots according to the order of the peaks here. 
%
%       comp: Other sessions you want to look at. 
%

%%  Preliminary. 
    %Load neuron map. 
    cd(mapMD.Location); load(fullfile(pwd,'batch_session_map.mat')); 
    
    %Compile session data. 
    sessions = [base,comp];
    dates = {sessions.Date};
    sessionNums = [sessions.Session];
    nSessions = length(sessions); 
    
    %Sort dates. 
    d = datenum(dates,'mm_dd_yyyy'); 
    sorted = sort(d); 
    [~,dateOrder] = ismember(d,sorted); 
    
    %Strings for titling. 
    dateTitles = dates; 
    for i=1:nSessions
        dateTitles{i}(3:3:6) = '-';
    end 
    
%% Find place cells. 
    cd(base.Location); load(fullfile(pwd,'PlaceMaps.mat'),'pval'); 
    PCs = find(pval>0.95); 

%% Find relevant indices in batch_session_map. 
    %Find the columns in MAP .
    [MAP,MAPcols] = FilterMAPDates(batch_session_map,dates,sessionNums);
    
    %Find row indices from the batch_session_map that correspond to time
    %cells in base. 
    MAProws = find(ismember(MAP(:,MAPcols(1)),PCs));
    
%% Plot. 
    f = figure('Position',[170 260 1280 460]); 
    for s=1:nSessions
        if s==1
            %Get the place cells. 
            neurons = MAP(MAProws,MAPcols(s)); 

            %Linearize trajectory and make matrix. 
            normRates = LinearizedPFs_treadmill(sessions(s));
            
            %Trim matrix to only include place cells. 
            normRates = normRates(neurons,:);
            
            %Find the peak responses and sort. 
            [~,inds] = max(normRates,[],2); 
            [~,order] = sort(inds); 
            plotme = normRates(order,:);
            
            %Plot. 
            subplot(1,nSessions,dateOrder(s));
            imagesc(plotme); colormap hot; 
            xlabel('Distance'); title(dateTitles{s});
            
        else
            %Get the same cells as in 'base'. 
            neurons = MAP(MAProws(order),MAPcols(s));      
            
            %Linearize trajectory and find spatial responses. 
            normRates = LinearizedPFs_treadmill(sessions(s));
               
            %Deal linearized place maps. 
            plotme = zeros(size(plotme));             
            for n=1:length(neurons)
                if neurons(n)
                    plotme(n,:) = normRates(neurons(n),:);         
                end
            end
            
            %Plot. 
            subplot(1,nSessions,dateOrder(s));
            imagesc(plotme); colormap hot; 
            xlabel('Distance'); title(dateTitles{s});
        end
        
        %Label. 
        if dateOrder(s)==1, ylabel('Neurons'); end
    end
    
    %For saving nicely as pdf.
    set(f,'PaperPositionMode','auto');         
    set(f,'PaperOrientation','landscape');
end