function TempSpatRasters(md,neurons)
%
%
%

%%
    path = md.Location;
    cd(path); 
    
    load('TimeCells.mat');
    nNeurons = length(neurons);
    delays = TodayTreadmillLog.delaysetting;
    complete = logical(TodayTreadmillLog.complete);
    nBins = unique(sum(~isnan(ratebylap(delays==T & complete,:,1)),2));
    bins = [1:.001:nBins]';
    t = linspace(0,T,length(bins));
     
    for n=1:nNeurons
        smoothfit = fit([1:nBins]',curves.tuning{neurons(n)}',...
            'smoothingspline');
        curves.smoothed{neurons(n)} = feval(smoothfit,bins); 
    end

    thisNeuron = 1; 
    keepgoing = true;
    sf = 0.1; 
    
    while keepgoing
        %Treadmill significance.
        [SIGX,SIGY] = significance(t,curves.sig{neurons(thisNeuron)},...
            curves.smoothed{neurons(thisNeuron)},bins);
        
        %Place field raster. 
        [placeRaster,placeCurve] = LinearizedPF_raster(md,'plotit',false,...
            'neurons',neurons(thisNeuron));
        
        f = figure(51);
        
        %Treadmill trial raster. 
        subplot(2,2,1);
        imagesc(0:T,1:sum(complete),ratebylap(:,:,neurons(thisNeuron)));
            colormap gray; freezeColors; 
            ylabel('Laps','fontsize',15); 
            title(['Neuron #',num2str(neurons(thisNeuron))]);
            set(gca,'fontsize',12,'ytick',[1,sum(complete)],'xtick',[]); 
            
        %Treadmill tuning curve. 
        subplot(2,2,3);
        plot(t,curves.smoothed{neurons(thisNeuron)},'color',[0 .5 .5],...
            'linewidth',5);
        hold on;
        Ylim = get(gca,'ylim');
        plot(SIGX,SIGY+Ylim(2)*sf,'ro','linewidth',4);
            xlim([0,T]);
            ylabel('Rate','fontsize',15);
            axis tight; 
            yLims = get(gca,'ylim');
            ylim([0,yLims(2)*(1+sf)]);
            set(gca,'tickdir','out','linewidth',4,'fontsize',12,'xtick',[0 5 10]);
            xlabel('Time (s)','fontsize',15);
            
        %Track trial raster.
        subplot(2,2,2); 
        imagesc(placeRaster);
        colormap hot;
        caxis([0 1.2]);
        set(gca,'xtick',[],'ytick',[1,sum(complete)],'linewidth',4,...
            'fontsize',12);
        
        %Track tuning curve. 
        subplot(2,2,4);
        plot(placeCurve,'color',[.58 .44 .86],'linewidth',5);
        xlabel('Linearized Distance (cm)','fontsize',15);
        axis tight; 
        yLims = get(gca,'ylim');
        ylim([0 yLims(2)]);
        set(gca,'linewidth',4','tickdir','out','fontsize',12);
           
        [keepgoing,thisNeuron] = scroll(thisNeuron,nNeurons,f); 
        close all;
    end
end

function [SIGX,SIGY] = significance(t,sigcurve,smoothedcurve,bins)
    %Get indices of significance.
    sigT = find(sigcurve);
    
    %Find the corresponding index in the smoothed bins.
    [~,inds] = ismember(sigT,bins);
    
    %Values to plot.
    SIGX = t(inds);
    SIGY = smoothedcurve(inds);
end