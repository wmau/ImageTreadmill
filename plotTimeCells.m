function plotTimeCells(animal,date,session,T)
%plotTimecells(animal,date,session,T)
%   
%   Plots single neuron responses in time during treadmill run. First
%   panel: spike-position trajectory. Second panel: Raster. Third panel:
%   Tuning curve. Press left/right arrows to scroll. Esc to exit. 
%
%   INPUTS
%       animal: Name of mouse (e.g., 'GCamp6f_45_treadmill').
%
%       date: Date of recording (e.g., '11_20_2015').
%
%       session: Session number. 
%   
%       T: Length of treadmill run. 
%

%%
    ChangeDirectory(animal,date,session);
        
    %Get treadmill data log. 
    TodayTreadmillLog = getTodayTreadmillLog(animal,date,session);
    TodayTreadmillLog = AlignTreadmilltoTracking(TodayTreadmillLog,TodayTreadmillLog.RecordStartTime);

    try
        load(fullfile(pwd,'TimeCells.mat')); 
    catch
        [TimeCells,ratebylap,curves,delays,x,y,time_interp,FT,T] = FindTimeCells(animal,date,session,T); 
    end
    
    [nNeurons,nFrames] = size(FT); 
    FT = logical(FT); 
    nBins = unique(sum(~isnan(ratebylap(delays==T,:,1)),2));
    
    %Get indices for treadmill runs. 
    inds = getTreadmillEpochs(TodayTreadmillLog,time_interp);
    temp = [];
    for thisLap=1:size(inds,1)
        if TodayTreadmillLog.complete(thisLap) && TodayTreadmillLog.delaysetting(thisLap) == T
            temp = [temp,inds(thisLap,1):inds(thisLap,2)];
        end
    end
    treadmillruns = logical(zeros(1,nFrames));
    treadmillruns(temp) = 1; 

%% Plot. 
    thisNeuron = 1;
    keepgoing = 1;
    sf = 0.1;
    bins = [1:0.001:nBins]';
    t = linspace(0,T,length(bins));     %Time vector, smoothed.
    tCI = linspace(0,T,length(curves.ci{TimeCells(thisNeuron)}(1,:)));  %Time vector for CI interpolation.
    
    %Simplify the matrix and get rid of nans if any. 
    ratebylap = ratebylap(delays==T,:,:);
    [~,c] = find(isnan(ratebylap(:,:,1)));
    ratebylap(:,c,:) = [];
    
    while keepgoing
        %Smooth the tuning curve.  
        smoothfit = fit([1:nBins]',curves.tuning{TimeCells(thisNeuron)}','smoothingspline');
        curves.smoothed{TimeCells(thisNeuron)} = feval(smoothfit,bins);
        
        %Get confidence intervals interpolated.
        shuffmean = mean(curves.shuffle{TimeCells(thisNeuron)});
        CImean = interp1(tCI,shuffmean,t,'pchip');
        CIu = interp1(tCI,curves.ci{TimeCells(thisNeuron)}(1,:),t,'pchip');
        CIl = interp1(tCI,curves.ci{TimeCells(thisNeuron)}(2,:),t,'pchip');
        
        %Get time bins of significance.
        [SIGX,SIGY] = significance(t,curves.sig{TimeCells(thisNeuron)},...
            curves.smoothed{TimeCells(thisNeuron)},bins);
        
        figure(50); 
        subplot(2,2,1);
            plot(x,y,x(treadmillruns & FT(TimeCells(thisNeuron),:)),y(treadmillruns & FT(TimeCells(thisNeuron),:)),'r.','MarkerSize',16);
            axis off; title(['Neuron #',num2str(TimeCells(thisNeuron))]);
        subplot(2,2,2);
            imagesc([0:T],[1:5:sum(delays==T)],ratebylap(:,:,TimeCells(thisNeuron)));
                colormap gray; ylabel('Laps');
        subplot(2,2,3:4);
            plot(t,curves.smoothed{TimeCells(thisNeuron)},'-r','linewidth',2);
            hold on; 
            plot(t,CImean,'-b','linewidth',2);
            plot(t,CIu,'--b',t,CIl,'--b');
            Ylim = get(gca,'ylim');
            plot(SIGX,SIGY+Ylim(2)*sf,'r*');
            hold off;
                xlabel('Time [s]'); ylabel('Rate'); 
                yLims = get(gca,'ylim');
                ylim([0, yLims(2)]);
                set(gca,'ticklength',[0 0]);
            
        figure(50);
            [~,~,key] = ginput(1); 

            if key == 29 && thisNeuron < length(TimeCells)
                thisNeuron = thisNeuron + 1; 
            elseif key == 28 && thisNeuron ~= 1
                thisNeuron = thisNeuron - 1; 
            elseif key == 27
                keepgoing = 0; 
                close(figure(50)); 
            end
            
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