function [AllI,nNeurons] = COMPILETIs(mds,celltype,infotype)
%[AllI,nNeurons] = COMPILETIs(mds,celltype,infotype)
%
%   

%% Set up.
    animals = unique({mds.Animal});
    nAnimals = length(animals); 

%% Compile
    AllI.stable = cell(1,nAnimals);
    AllI.unstable = cell(1,nAnimals);
    nNeurons.stable = zeros(1,nAnimals); 
    nNeurons.unstable = zeros(1,nAnimals); 
    for a = 1:nAnimals
        nStable = 0; 
        nUnstable = 0;
        AllI.stable{a} = [];
        AllI.unstable{a} = [];
        
        %Get all the sessions for this animal.
        ssns = find(strcmp(animals{a},{mds.Animal})); 
        for s = ssns(1:end-1)
            cd(mds(s).Location);
            
            if strcmp(infotype,'time');
                load('TemporalInfo.mat','I');
            elseif strcmp(infotype,'place')
                load('SpatialInfo.mat','spatialI');
                I = spatialI';
            end
            load('TimeCells.mat','TimeCells');
            load('PlaceMaps.mat','pval');
            load('PFstats.mat','PFnumhits','MaxPF');
            
            %Get all time cells with a viable place field. 
            idx = sub2ind(size(PFnumhits), 1:size(PFnumhits,1), MaxPF);
            noi = intersect(TimeCells,find(pval>.95 & PFnumhits(idx) > 4));
            
            if strcmp(celltype,'time')
                %Get correlation coefficients and p-values. 
                corrStats = CorrTrdmllTrace(mds(s),mds(s+1),noi);
                tuningStatus = TCRemap(mds(s),mds(s+1));

                %Stable time cells based on correlation and non-shifting time
                %field.
                stable = find(corrStats(:,2) < .01 & tuningStatus(:,2) > 0);
                unstable = find(corrStats(:,2) > .01 | tuningStatus(:,2) < 1);
         
            elseif strcmp(celltype,'place')
                %Get the correlation coefficients and p-values.
                corrStats = CorrPlaceFields(mds(s),mds(s+1),noi);

                %Stable place cells based on correlation p-value.
                stable = find(corrStats(:,2) <= .01); 
                unstable = find(corrStats(:,2) > .01);     
            end
           
            %Get number of stable and unstable cells. Add this one per
            %session but track the number per animal. 
            nStable = nStable + length(stable); 
            nUnstable = nUnstable + length(unstable); 
            
            %Get the temporal information values. 
            AllI.stable{a} = [AllI.stable{a}; I(stable)];
            AllI.unstable{a} = [AllI.unstable{a}; I(unstable)];
        end
        
        %Get number of neurons that are (un)stable per animal. 
        nNeurons.stable(a) = nStable; 
        nNeurons.unstable(a) = nUnstable;
    end
    
end
        
        