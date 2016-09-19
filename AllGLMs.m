function [neurons,all,nocells] = AllGLMs(md)
%[neurons,all,nocells,automated] = AllGLMs(md)
%
%

%%
    cd(md.Location);
    load('graphData_p.mat','A');
    [nNeurons,neurons] = nNeuronsActiveonTM(md);

    [X,y] = SourceSinkGLMSetUp(md,neurons,1); %Don't delete extra output here. 
    
    nn=1;
    all = cell(1,nNeurons);
    notime = cell(1,nNeurons);
    nocells = cell(1,nNeurons);
    automated = cell(1,nNeurons);
    p = ProgressBar(nNeurons);
    for n=neurons
        [~,y] = SourceSinkGLMSetUp(md,[],n);
        fitTbl = X;                     %Get design matrix.
        fitTbl(:,n+1) = [];             %Take out the row corresponding to identity.
        fitTbl(:,end+1) = table(y);     %Add in response variable. 
        
        %Do the fit. 
        [all{nn},notime{nn},nocells{nn}] = SourceSinkGLM(fitTbl);
        
        nn=nn+1;
        
        p.progress;
    end
    p.stop;
end