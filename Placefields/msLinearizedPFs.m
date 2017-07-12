function msLinearizedPFs(base,comp,varargin)
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
    p = inputParser;
    p.addRequired('base',@(x) isstruct(x));
    p.addRequired('comp',@(x) isstruct(x));
    p.addParameter('noTreadmill',true,@(x) islogical(x));
    
    p.parse(base,comp,varargin{:});
    
    noTreadmill = p.Results.noTreadmill;
    
    %Compile session data. 
    sessions = [base,comp];
    dates = {sessions.Date};
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
    PCcrit = .01;
    PCs = getPlaceCells(base,PCcrit);

%% Find relevant indices in batch_session_map. 
    matches = msMatchCells(mapMD,sessions,PCs,false);
    nNeurons = size(matches,1);
    
%% Plot. 
    f = figure('Position',[170 260 260*nSessions 460]); 
    for s=1:nSessions
        neurons = matches(:,s);
        nNeurons = length(neurons);
        
        if s==1
            %Get the place cells. 

            %Linearize trajectory and make matrix. 
            [~,normRates,~,~,~,edges] = LinearizedPFs_treadmill(sessions(s),'plotit',false);
            
            %Trim matrix to only include place cells. 
            normRates = normRates(neurons,:);
            
            %Find the peak responses and sort. 
            [~,inds] = max(normRates,[],2); 
            [~,order] = sort(inds); 
            plotme = normRates(order,:);
            
            [bad,~] = find(isnan(plotme));
            plotme(bad,:) = 0;
            
            %Plot. 
            subplot(1,nSessions,dateOrder(s)); 
            imagesc(plotme); colormap hot; hold on;
            plot(inds(order),1:nNeurons,'color',[.58 .44 .86],'linewidth',5);
            set(gca,'ytick',[1,nNeurons]); axis tight;
            if noTreadmill
                xlim([22 80]);
                set(gca,'xtick',[22 80],'xticklabel',[0 140]);
            end
            xlabel('Track position (cm)'); title(dateTitles{s}); 
            
        else
            neurons = neurons(order);
            %Linearize trajectory and find spatial responses. 
            [~,normRates] = LinearizedPFs_treadmill(sessions(s),'plotit',false);
               
            %Deal linearized place maps. 
            plotme = zeros(nNeurons,length(edges)-1);             
            for n=1:length(neurons)
                if ~isnan(neurons(n)) && neurons(n)>0
                    plotme(n,:) = normRates(neurons(n),:);
                else, plotme(n,:) = 0;
                end               
            end
            [bad,~] = find(isnan(plotme));
            plotme(bad,:) = 0;
            
            %Plot. 
            subplot(1,nSessions,dateOrder(s)); 
            imagesc(plotme); colormap hot; hold on;
            plot(inds(order),1:nNeurons,'color',[.58 .44 .86],'linewidth',5,...
                'linestyle',':');
            set(gca,'ytick',[1,nNeurons]); axis tight;
            if noTreadmill
                xlim([22 80]);
                set(gca,'xtick',[22 80],'xticklabel',[0 140]);
            end
            xlabel('Distance'); title(dateTitles{s});
        end
        
        %Label. 
        if dateOrder(s)==1, ylabel('Neurons'); end
    end
    
    %For saving nicely as pdf.
    set(f,'PaperPositionMode','auto');         
    set(f,'PaperOrientation','landscape');
end