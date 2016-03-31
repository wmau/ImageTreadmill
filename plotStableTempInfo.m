function [sigI,nsigI] = plotStableTempInfo(mapMD,MD1,MD2) 
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
    [~,tCorr] = PlaceTimeCorr(mapMD,MD1,MD2,TimeCells); 

%% Get stable neurons. 
    load(fullfile(MD1.Location,'TemporalInfo.mat'),'I'); 
    
    sig = tCorr(:,2) < 0.05;
    nsig = tCorr(:,2) > 0.05;
    
    sigI = I(sig);
    nsigI = I(nsig);

%% Plot.
    %Histogram. 
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