function [I,sigI,nsigI,tCorr,pCorr] = plotStableTempInfo(mapMD,MD1,MD2,type,plotit) 
%plotStableTempInfo(mapMD,MD1,MD2) 
%
%   Investigates the temporal information content of neurons based on their
%   stability the next day. Answers the question: Do time cells that retain
%   their tuning contain more temporal information? 
%
%   INPUTS
%       mapMD: MD entry where batch_session_map lives.
%
%       MD1: Looks at temporal information on this session. 
%
%       MD2: Look at whether time cells in MD1 are stable on this session.
%       Stable neurons are those that have a significant correlation
%       between tuning curves. 
%
%   OUTPUT
%       Histogram and ECDF of temporal information grouped by stable and
%       non-stable neurons. 
%

%% Do correlations.
    load(fullfile(MD1.Location,'TimeCells.mat'),'TimeCells'); 
    [pCorr,tCorr] = PlaceTimeCorr(mapMD,MD1,MD2,TimeCells); 
    
%% Alternatively, look for shifts in tuning centroids. 
    s = TimeCellRemapRate(mapMD,MD1,MD2,[10 10]);
    sig = false(length(pCorr),1);
    nsig = false(length(pCorr),1);
    
%% Get stable neurons. 
    load(fullfile(MD1.Location,'TemporalInfo.mat'),'I'); 
    load(fullfile(MD1.Location,'PlaceMaps.mat'),'pval'); 
    
    if strcmp(type,'time')
        sig = tCorr(:,2) < 0.05;
        nsig = tCorr(:,2) > 0.05;
        
        %sig(TimeCells(s(:,2)==1)) = true;
        %nsig(TimeCells(s(:,2)~=1)) = true; 
    elseif strcmp(type,'place')
        sig = pCorr(:,2) < 0.05;
        nsig = pCorr(:,2) > 0.05;
    end
    
    sigI = I(sig);
    nsigI = I(nsig);

%% Plot.
    %Histogram. 
    if plotit
        figure('position',[390 320 800 385]);
        subplot(1,2,1); hold on;
            %Get edges. 
            [~,edges] = histcounts(I(sig),linspace(0,max(I(sig | nsig)),25));
            [n,b] = hist(sigI,edges);     %Bin.
            n = n./length(sigI);          %Normalize. 
            stairs(b,n);
            [n,b] = hist(nsigI,edges);    %Bin.
            n = n./length(nsigI);         %Normalize. 
            stairs(b,n); 
                xlim([0 max(I(sig | nsig))]);
                set(gca,'ticklength',[0 0]);
                xlabel('Temporal Information [bits/s]');
                ylabel('Proportion'); 
                hold off;

        %ECDF. 
        subplot(1,2,2); hold on;
            ecdf(sigI); 
            ecdf(nsigI); 
                xlim([0 max(I(sig | nsig))]);
                legend('Stable','Unstable','location','southeast'); 
                xlabel('Temporal Information [bits/s]'); 
                ylabel('Proportion'); 
                set(gca,'ticklength',[0 0]);
                hold off;
    end
end