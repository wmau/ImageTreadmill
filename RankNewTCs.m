function [S1Peaks,S2Peaks] = RankNewTCs(base,comp,varargin)
%
%
%

%% 
    p = inputParser; 
    p.addRequired('base',@(x) isstruct(base));
    p.addRequired('comp',@(x) isstruct(comp)); 
    p.addParameter('plotit',true,@(x) islogical(x));
    
    p.parse(base,comp,varargin{:});
    
    plotit = p.Results.plotit;
    
%% 
    cd(base.Location);
    load('TimeCells.mat','T');
    [S1,S2] = getNewTimeCells(base,comp); 
    
    [~,S2order,peaks] = PastalkovaPlot(comp,'TimeCells',S2,'plotit',false); 
    peaks = peaks/4;
    [~,S1order] = PastalkovaPlot(base,'TimeCells',S1,'plotit',false);
    PastalkovaPlot(base,'TimeCells',S1,'order',S2order,'plotit',false); 
    
    if plotit
        figure('Position',[-1250 250 560 420]);
        subplot(1,2,2);
            PastalkovaPlot(comp,'TimeCells',S2);
            title('Day 2');

            nTCs = length(S2order); 
            hold on;
            plot(peaks,[1:nTCs],'r','linewidth',3);
        subplot(1,2,1); 
            PastalkovaPlot(base,'TimeCells',S1);
            title('Day 1');

            hold on;
            plot(peaks,[1:nTCs],'r','linewidth',3);

        figure;
        scatter(S1order,S2order,'.');
        [R,p] = corr(S1order,S2order,'type','spearman');
        title(['R = ',num2str(R), ' p = ',num2str(p)]);
    end
    
    [~,temp] = getTimePeak(base); 
    S1Peaks = temp(S1);
    
    [~,temp] = getTimePeak(comp);
    S2Peaks = temp(S2); 
end