function [rawtrdmll,difftrdmll,tracetrdmll] = TreadmillTraces(md)
%[rawtrdmll,difftrdmll,tracetrdmll] = TreadmillTraces(md)
%
%

%%
    cd(md.Location);
    load('TimeCells.mat','TimeCells','TodayTreadmillLog','T');
    load('Pos_align.mat','rawtrace','difftrace','trace');
    
    [inds,nRuns] = TrimTrdmllInds(TodayTreadmillLog,T); 
    nNeurons = size(rawtrace,1);
    
    rawtrdmll = nan(nRuns,T*20,nNeurons);
    difftrdmll = nan(nRuns,T*20,nNeurons);
    tracetrdmll = nan(nRuns,T*20,nNeurons);
    window = ones(10,1)/10;
    for n=1:nNeurons
        rawtrdmll(:,:,n) = buildRasterTrace(inds,rawtrace,n);
        difftrdmll(:,:,n) = buildRasterTrace(inds,difftrace,n);
        for l=1:nRuns
            difftrdmll(l,:,n) = convtrim(difftrdmll(l,:,n),window);
        end
        tracetrdmll(:,:,n) = buildRasterTrace(inds,trace,n);
    end
    
    save('TreadmillTraces','rawtrdmll','difftrdmll','tracetrdmll');
end