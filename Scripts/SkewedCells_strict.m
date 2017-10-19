clear;
loadMD;

fulldataset = [MD(292:303) MD(305:308)];            %Sessions.
nSessions = length(fulldataset);                    %Number of sessions.
cellType = 'timecells';                             %Only time cells work right now.
rasterType = 'time';                                %Look at temporal responses. 
B = 500;                                             %Number of trial shuffles. 

[skew,sig,strictSig,even,odd] = deal(cell(nSessions,1)); 
for s=1:nSessions
    disp(['Analyzing ',fulldataset(s).Animal,' ',fulldataset(s).Date]);
    
    %Compute skewnesses for all trials and even/odd trials. 
    [skew{s},even{s},odd{s}] = getAllSkewnesses(fulldataset(s),'cellType',...
        cellType,'rasterType',rasterType,'subsample',true); 
    
    %Shuffle trials and see what skewness scores you get. 
    shuffled = cell(B,1); 
    for i=1:B
        shuffled{i} = getAllSkewnesses(fulldataset(s),'cellType',cellType,...
            'rasterType',rasterType,'shuffle',true); 
    end
    
    %Get time cells and then check how much their skewnesses exceed that
    %given chance. 
    neurons = AcquireTimePlaceCells(fulldataset(s),cellType)'; 
    [sig{s},strictSig{s}] = deal(nan(size(skew{s}))); 
    for n=neurons
        %Get the distribution of skewnesses from trial-shuffled rasters. 
        distribution = cellfun(@(x) x(n),shuffled); 
        
        %P-value, the percentage of skewnesses that exceed chance. 
        p = sum(skew{s}(n) > distribution)/B; 
        even_p = sum(even{s}(n) > distribution)/B;  %Only even trials.
        odd_p = sum(odd{s}(n) > distribution)/B;    %Only odd trials. 
       
        %Mark cells as significantly trial-skewed. 
        if p < 0.025            %If p < 0.025, they are early-skewed cells. 
            sig{s}(n) = 1; 
            
            %Strict criterion where both even and odd trials must be
            %significantly different from chance. 
            if even_p < 0.025 && odd_p < 0.025
                strictSig{s}(n) = 1; 
            else
                strictSig{s}(n) = 0;
            end
            
        elseif p > 0.975        %If p > 0.975, they are late-skewed cells. 
            sig{s}(n) = 2; 
            
            if even_p > 0.975 && odd_p > 0.975
                strictSig{s}(n) = 2; 
            else
                strictSig{s}(n) = 0;
            end
        end
    end
end

nEarly = sum(cellfun(@(x) sum(x==1),strictSig));
nLate = sum(cellfun(@(x) sum(x==2),strictSig));

