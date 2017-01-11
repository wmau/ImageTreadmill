function msLinearizedPFs(base,comp)
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
    
    %Compile session data. 
    sessions = [base,comp];
    dates = {sessions.Date};
    sessionNums = [sessions.Session];
    nSessions = length(sessions); 
    
    mapMD = getMapMD(sessions);
    
    %Load neuron map. 
    cd(mapMD.Location); load(fullfile(pwd,'batch_session_map.mat')); 
    
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
    cd(base.Location); 
    load(fullfile(pwd,'Placefields.mat'),'pval'); 
    load(fullfile(pwd,'PlacefieldStats.mat'),'PFnHits','bestPF');
    load(fullfile(pwd,'SpatialInfo.mat'),'MI');
    idx = sub2ind(size(PFnHits), 1:size(PFnHits,1), bestPF');
    PCcrit = .01;
    
    PCs = find(pval<PCcrit & MI'>0 & PFnHits(idx)>4);

%% Find relevant indices in batch_session_map. 
    matches = msMatchCells(mapMD,sessions,PCs,false);
    nNeurons = size(matches,1);
    
%% Plot. 
    f = figure('Position',[170 260 1280 460]); 
    for s=1:nSessions
        neurons = matches(:,s);
        
        if s==1
            %Get the place cells. 

            %Linearize trajectory and make matrix. 
            [~,normRates,~,~,~,edges] = LinearizedPFs_treadmill(sessions(s));
            
            %Trim matrix to only include place cells. 
            normRates = normRates(neurons,:);
            
            %Find the peak responses and sort. 
            [~,inds] = max(normRates,[],2); 
            [~,order] = sort(inds); 
            plotme = normRates(order,:);
            
            %Plot. 
            subplot(1,nSessions,dateOrder(s));
            imagesc(edges,1:nNeurons:nNeurons,plotme); colormap hot; 
            xlabel('Distance'); title(dateTitles{s});
            
        else
            neurons = neurons(order);
            %Linearize trajectory and find spatial responses. 
            [~,normRates] = LinearizedPFs_treadmill(sessions(s));
               
            %Deal linearized place maps. 
            plotme = zeros(size(plotme));             
            for n=1:length(neurons)
                if ~isnan(neurons(n)) && neurons(n)>0
                    plotme(n,:) = normRates(neurons(n),:);
                end               
            end
            
            %Plot. 
            subplot(1,nSessions,dateOrder(s));
            imagesc(edges,1:nNeurons:nNeurons,plotme); colormap hot; 
            xlabel('Distance'); title(dateTitles{s});
        end
        
        %Label. 
        if dateOrder(s)==1, ylabel('Neurons'); end
    end
    
    %For saving nicely as pdf.
    set(f,'PaperPositionMode','auto');         
    set(f,'PaperOrientation','landscape');
end