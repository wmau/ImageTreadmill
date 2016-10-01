function [neurons,nocells,automated] = SourceOnlyGLMs(md,lags,tracetype)
%[neurons,all,nocells,automated] = AllGLMs(md)
%
%

%%
    cd(md.Location);
    load('XCorr.mat','A');
    [nNeurons,neurons] = nNeuronsActiveonTM(md);

    %Set up design matrix. 
    [X,y] = SourceSinkGLMSetUp(md,neurons,1,lags,tracetype); %Don't delete extra output here.
    
    nn=1;
    all = cell(1,nNeurons);
    notime = cell(1,nNeurons);
    nocells = cell(1,nNeurons);
    automated = cell(1,nNeurons);
    p = ProgressBar(nNeurons);
    for n=neurons
        [~,y] = SourceSinkGLMSetUp(md,[],n,0,tracetype);
        
        sources = find(A(:,n)); cols = 1;
        for s=sources'
            cols = [cols, find(~cellfun('isempty',regexp(X.Properties.VariableNames,['n',num2str(s),'lag'])))];
        end
        
        if ~isempty(cols)
            fitTbl = X(:,cols);                 %Get design matrix.  
            %fitTbl = X;
            %fitTbl(:,['n',num2str(n)]) = [];    %Get rid of identity.        
            fitTbl(:,end+1) = table(y);         %Add in response variable. 

            %Do the fit. 
            [nocells{nn},automated{nn}] = SourceSinkGLM(fitTbl,tracetype);
        end
        
        nn=nn+1;
        
        p.progress;
    end
    p.stop;
end