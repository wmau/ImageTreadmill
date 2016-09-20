function [neurons,all,notime,nocells] = AllGLMs(md,tracetype)
%[neurons,all,nocells,automated] = AllGLMs(md)
%
%

%%
    cd(md.Location);
    load('graphData_p.mat','A');
    [nNeurons,neurons] = nNeuronsActiveonTM(md);

    [X,y] = SourceSinkGLMSetUp(md,neurons,1,tracetype); %Don't delete extra output here.
    
    nn=1;
    all = cell(1,nNeurons);
    notime = cell(1,nNeurons);
    nocells = cell(1,nNeurons);
    automated = cell(1,nNeurons);
    p = ProgressBar(nNeurons);
    for n=neurons
        [~,y] = SourceSinkGLMSetUp(md,[],n,tracetype);
        
        sources = find(A(:,n));
        [~,Xind] = ismember(sources,neurons);   
        
        %fitTbl = X(:,[1;Xind+1]);           %Get design matrix.  
        fitTbl = X;
        fitTbl(:,['n',num2str(n)]) = [];     %Get rid of identity.
        fitTbl(:,end+1) = table(y);         %Add in response variable. 
        
        %Do the fit. 
        [all{nn},notime{nn},nocells{nn}] = SourceSinkGLM(fitTbl,tracetype);
        
        nn=nn+1;
        
        p.progress;
    end
    p.stop;
end