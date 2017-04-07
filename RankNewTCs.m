function [S1Peaks,S2Peaks,S1] = RankNewTCs(base,comp,varargin)
%[S1Peaks,S2Peaks] = RankNewTCs(base,comp,varargin)
%
%   Gets peaks of new time cells and matches them to the peaks of
%   themselves the previous day. 
%
%   INPUTS
%       base: Day before. 
%
%       comp: Day neurons became time cells. 
%
%   varargin: 
%       plotit: Logical, whether or not to Pastalkova plot and scatter
%       peaks and ranks.
%

%% Parse inputs. 
    p = inputParser; 
    p.addRequired('base',@(x) isstruct(base));
    p.addRequired('comp',@(x) isstruct(comp)); 
    p.addParameter('plotit',true,@(x) islogical(x));
    
    p.parse(base,comp,varargin{:});
    
    plotit = p.Results.plotit;
    teal = [0 .5 .5];
    
%% Get new time cells. 
    cd(base.Location);
    load('TimeCells.mat','T');
    [S1,S2] = getNewTimeCells(base,comp); 
    
%% Meat of the function.
    %Order the new time cells. 
    [~,S2order,peaks] = PastalkovaPlot(comp,'TimeCells',S2,'plotit',false); 
    
    %Conversion from bin size. 
    peaks = peaks/4;
    
    %Get the order of the new time cells the day before. 
    [~,S1order] = PastalkovaPlot(base,'TimeCells',S1,'plotit',false);
    
    %Get peaks of the new time cells from session 1 in the same order as session 2.
    [~,temp] = getTimePeak(base); 
    S1Peaks = temp(S1); 
    [~,temp] = getTimePeak(comp);
    S2Peaks = temp(S2); 

%% Plotting.    
    if plotit
        figure('Position',[-1250 250 560 420]);
        %Sorted day 2 new time cells. 
        subplot(1,2,2); hold on;
            [~,~,~,TCDay] = PastalkovaPlot(comp,'TimeCells',S2); 
            nTCs = length(S2order); 
            plot(peaks,[1:nTCs],'color',teal,'linewidth',3);

            title('Day 0');
            set(gca,'ytick',[1,nTCs]);
        
            TCMaxes = max(TCDay,[],2); 
        %New time cells on day 1 sorted in the same order as day 2. 
        subplot(1,2,1); hold on;
            [~,~,~,preTCDay] = PastalkovaPlot(base,'TimeCells',S1,'order',S2order,'plotit',false);
            T = size(preTCDay,2)/4; 
            imagesc([0:T],[1:nTCs],preTCDay./TCMaxes);
            set(gca,'ydir','reverse'); axis tight;
            plot(peaks,[1:nTCs],'color',teal,'linewidth',3,'linestyle','--');
            
            set(gca,'ytick',[1,nTCs]);
            title('Day -1');  
            ylabel('Neurons');
       
        %Scatter plot of peaks. 
%         figure;
%         scatter(S1Peaks,S2Peaks,'.');
%         [R,p] = corr(S1Peaks,S2Peaks,'type','spearman');
%         title(['R = ',num2str(R), ' p = ',num2str(p)]);
%         xlabel('Day 1 Peak'); 
%         ylabel('Day 2 Peak');
    end
    
    
end