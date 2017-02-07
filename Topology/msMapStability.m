function msMapStability(base,comp,stabilityType)
%
%
%

%%
    cd(base.Location); 
    load('FinalOutput.mat','NumNeurons'); 
    PCcrit = .01;
    switch stabilityType
        case 'time'
            neurons = getTimeCells(base);
            corrs = CorrTrdmllTrace(base,comp,neurons);
            color = [0 .5 .5];
        case 'place'
            neurons = getPlaceCells(base,PCcrit); 
            corrs = CorrPlaceFields(base,comp,neurons);
            color = [.58 .44 .86];
    end

    stblcrit = .01/length(neurons);
    
    stable = intersect(find(corrs(:,2) < stblcrit),neurons); 
    unstable = intersect(find(corrs(:,2) > stblcrit | isnan(corrs(:,2))),neurons); 
    
    PlotNeurons(base,1:NumNeurons,[0 0 0 .2],1);
    hold on; 
    PlotNeurons(base,stable,color,2);
    PlotNeurons(base,unstable,[1 0 0 .5],2);

            
    end