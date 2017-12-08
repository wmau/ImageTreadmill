c = 1;
[agreeStable,agreeUnstable,totalStable,totalUnstable] = deal(zeros(1,14));
for i=1:4
    nSessions = length(stableCorr{i});
    
    for s=1:nSessions
        agreeStable(c) = length(intersect(stableCorr{i}{s},stableShift{i}{s}));
        agreeUnstable(c) = length(intersect(unstableCorr{i}{s},unstableShift{i}{s}));
        
        disp(['Correlation unstable: ',num2str(length(unstableCorr{i}{s}))]);
        disp(['Shift unstable: ',num2str(length(unstableShift{i}{s}))]);
        
        totalStable(c) = length(stableShift{i}{s});
        totalUnstable(c) = length(unstableShift{i}{s}); 
        
        c = c+1;
    end
end

mean([agreeStable./totalStable agreeUnstable./totalUnstable])
std([agreeStable./totalStable agreeUnstable./totalUnstable])./sqrt(14)
(sum(agreeStable)+sum(agreeUnstable))/(sum(totalStable)+sum(totalUnstable))