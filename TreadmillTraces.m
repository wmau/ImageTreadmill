function [RawTrdmll,DFDTTrdmll,LPtrdmll] = TreadmillTraces(md)
%[rawtrdmll,difftrdmll,tracetrdmll] = TreadmillTraces(md)
%
%

%%
    cd(md.Location);
    load('TimeCells.mat','TimeCells','TodayTreadmillLog','T');
    load('Pos_align.mat','RawTrace','DFDTtrace','LPtrace');
    
    [inds,nRuns] = TrimTrdmllInds(TodayTreadmillLog,T); 
    nNeurons = size(RawTrace,1);
    
    RawTrdmll = nan(nRuns,T*20,nNeurons);
    DFDTTrdmll = nan(nRuns,T*20,nNeurons);
    LPtrdmll = nan(nRuns,T*20,nNeurons);
    window = ones(10,1)/10;
    for n=1:nNeurons
        RawTrdmll(:,:,n) = buildRasterTrace(inds,RawTrace,n);
        DFDTTrdmll(:,:,n) = buildRasterTrace(inds,DFDTtrace,n);
        for l=1:nRuns
            DFDTTrdmll(l,:,n) = convtrim(DFDTTrdmll(l,:,n),window);
        end
        LPtrdmll(:,:,n) = buildRasterTrace(inds,LPtrace,n);
    end
    
    save('TreadmillTraces','RawTrdmll','DFDTTrdmll','LPtrdmll');
end