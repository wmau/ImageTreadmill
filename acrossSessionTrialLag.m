function acrossSessionTrialLag(R,lapNum,sessionNum)
% acrossSessionTrialLag(R,lapNum,sessionNum)
%
% Look at lag effects across days AND sessions. e.g., trial 1 on day 1 vs
% trial 2 on day 2, trial 3 on day 2, trial 4 on day 2, etc.

%% 
    nSessions = max(sessionNum);
   
    [m,sem,diags] = deal(cell(nSessions));
    for s1=1:nSessions
        row = sessionNum==s1;
        
        for s2=s1:nSessions
            col = sessionNum==s2;
            
            %Get the block of correlation coefficients corresponding to
            %that day-to-day comparison. 
            block = R(row,col); 
            
            [m{s1,s2},sem{s1,s2},diags{s1,s2}] = collapseByLag(block);
        end
    end
   
%     figure;
%     M = m';
%     SEM = sem';
%     for i=1:25
%         subplot(5,5,i);
% 
%         if ~isempty(M{i})
%             errorbar(M{i},SEM{i});
%         end
%     end
    
    nElements = cellfun('length',m);
    sessionLag = cell(nSessions,1);
    figure; hold on;
    for s=1:nSessions
        row = 1:nSessions-s+1;
        col = s:nSessions;
        
        inds = sub2ind([nSessions,nSessions],row,col);
        
        trimLength = min(nElements(inds));
        
        sessionLag{s} = cell(trimLength,1);
        for l=1:trimLength          
            toInsert = cell2mat(cellfun(@(x) x(l), {diags{inds}}));
            sessionLag{s}{l} = [sessionLag{s}{l} toInsert];
        end
        
        a = cellfun(@nanmean,sessionLag{s});
        b = cellfun(@nanstd,sessionLag{s}) ./ sqrt(cellfun('length',sessionLag{s}));
       
        errorbar(a,b,'capsize',0);
    end
    
end