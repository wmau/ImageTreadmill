function [corrs,cellID] = PairwiseTrdmllCorrs(mds,neurons)
%PairwiseTrdmllCorrs(mds,neurons)
%
%   Computes the correlation coefficient between the tuning curves of
%   adjacent sessions specified in mds. 
%
%   INPUTS
%       mds: Session entries for a single animal. 
%
%       neurons: Neurons 

%%
    nSessions = length(mds)-1;
    nNeurons = length(neurons); 
    
    %Map the neurons.
    mapMD = getMapMD(mds);
    map = msMatchCells(mapMD,mds,neurons,false);
    cellID = map(:,1);
    
    %Preallocate. 
    corrs = cell(nNeurons,1);  
    for s=1:nSessions
        %Only take the neurons that mapped. 
        good = map(:,s) > 0;
        
        %Correlate the neurons that mapped. 
        c = CorrTrdmllTrace(mds(s),mds(s+1),map(good,s));
        
        %Compile.
        for n=1:nNeurons
            %If neuron got mapped...
            if map(n,s) > 0
                %Get its correlation coefficient. 
                corrs{n}(s,:) = c(map(n,s),:); 
            else, corrs{n}(s,:) = [0 1]; 
            end 
            
        end
    end
end