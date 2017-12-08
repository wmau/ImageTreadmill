function matrix = weirdassBinningMethod(R,binSize,nBins,slidingWindow)
%
%
%

%%    
    atTheEnd = false;
    nTrials = size(R,1);
    matrix = cell(nBins);
    i = 1;
    while ~atTheEnd        
        if i==1
            startInd = 0; 
            endInd = nBins*binSize; 
        else
            if slidingWindow 
                startInd = startInd + 1; 
                endInd = endInd + 1; 
            else
                startInd = endInd; 
                endInd = startInd + nBins*binSize;
            end
        end
        
        %Vector containing the limits on trial blocks. 
        trialBlocks = startInd:binSize:endInd; 
        
        if trialBlocks(end) > nTrials
            trialBlocks = startInd:binSize:nTrials;
            nBins = length(trialBlocks)-1;
            atTheEnd = true;
        end
            
        for s1=1:nBins
            row = trialBlocks(s1)+1:trialBlocks(s1+1);
            
            for s2=s1:nBins
                col = trialBlocks(s2)+1:trialBlocks(s2+1);

                toAppend = R(row,col);
                toAppend = toAppend(:);
                
                matrix{s1,s2} = [matrix{s1,s2}; toAppend]; 
                
                if s2~=s1
                    matrix{s2,s1} = [matrix{s2,s1}; toAppend]; 
                end
            end
        end
               
        i = i+1;
    end
    
    
end