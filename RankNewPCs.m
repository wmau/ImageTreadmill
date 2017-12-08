function [peaks1,peaks2,S1] = RankNewPCs(base,comp,varargin)
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
    purple = [.58 .44 .86];
    
%% Get new time cells. 
    cd(base.Location);
    [S1,S2] = getNewPlaceCells(base,comp); 
    
%% Meat of the function.
    %Order the new time cells. 
    [~,~,~,S2order,~,~,peaks1] = LinearizedPFs_treadmill(comp,'PlaceCells',S2,'plotit',false); 
    
    %Get the order of the new time cells the day before. 
    [~,~,~,S1order,~,~,peaks2] = LinearizedPFs_treadmill(base,'PlaceCells',S1,'plotit',false);

%% Plotting.    
    if plotit
        figure('Position',[-1250 250 560 420]); 
        %Sorted day 2 new time cells. 
        subplot(1,2,2); hold on;
            PCDay = LinearizedPFs_treadmill(comp,'PlaceCells',S2);
            nPCs = length(S2order);
            plot(peaks1,[1:nPCs],'color',purple,'linewidth',3);
            
            title('Day 0');
            set(gca,'ytick',[1,nPCs]);  
            set(gca,'xtick',[1,size(PCDay,2)]);
            
            PCMaxes = max(PCDay,[],2); 
        %New time cells on day 1 sorted in the same order as day 2. 
        subplot(1,2,1); hold on;
            prePCDay = LinearizedPFs_treadmill(base,'PlaceCells',S1,'order',S2order,'plotit',false);
            imagesc([1:size(prePCDay,2)],[1:nPCs],prePCDay./PCMaxes); axis tight;
            set(gca,'ydir','reverse'); axis tight;
            plot(peaks1,[1:nPCs],'color',purple,'linewidth',3,'linestyle','--');
            
            title('Day -1');
            set(gca,'ytick',[1,nPCs]);
            set(gca,'xtick',[1,size(PCDay,2)]);
            ylabel('Neurons');

        %Scatter plot of orders.  
%         figure; hold on;
%         scatter(peaks1,peaks2,'.');
%         [R,p] = corr(peaks1,peaks2,'type','spearman');
%         title(['R = ',num2str(R), ' p = ',num2str(p)]);
%         xlabel('Day 1 Rank'); 
%         ylabel('Day 2 Rank');
        
    end
    
    
end