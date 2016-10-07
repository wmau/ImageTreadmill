function [neurons,nocells,automated] = AllGLMs(md,lags,tracetype)
%[neurons,all,nocells,automated] = AllGLMs(md)
%
%

%%
    cd(md.Location);
    [nNeurons,neurons] = nNeuronsActiveonTM(md);

    %Set up design matrix. 
    [X,y] = SourceSinkGLMSetUp(md,neurons,1,lags,tracetype); %Don't delete extra output here.
    
    nn=1;
    nocells = cell(1,nNeurons);
    automated = cell(1,nNeurons);
    p = ProgressBar(nNeurons);
    for n=neurons
        [~,y] = SourceSinkGLMSetUp(md,[],n,0,tracetype);
                                 
        identity = find(~cellfun('isempty',regexp(X.Properties.VariableNames,['n',num2str(n),'lag'])));
        fitTbl = X;                         %Get design matrix.
        fitTbl(:,identity) = [];            %Get rid of identity.        
        fitTbl(:,end+1) = table(y);         %Add in response variable. 

        %Do the fit. 
        [nocells{nn},automated{nn}] = SourceSinkGLM(fitTbl,tracetype);

        nn=nn+1;
        
        p.progress;
    end
    p.stop;
end