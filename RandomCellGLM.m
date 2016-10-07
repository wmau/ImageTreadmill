function [neurons,randcells] = RandomCellGLM(md,lags,autoGLMs,tracetype)
%[neurons,all,nocells,automated] = AllGLMs(md)
%
%

%%
    cd(md.Location);
    load('XCorr.mat','A');
    [~,snks] = find(A); 
    snks = unique(snks)';
    [~,neurons] = nNeuronsActiveonTM(md);
    nNeurons = size(A,1);

    %Set up design matrix. 
    [X,y] = SourceSinkGLMSetUp(md,neurons,1,lags,tracetype); %Don't delete extra output here.
    
    nn=1;
    B = 500;
    randcells = cell(B,nNeurons);
    p = ProgressBar(length(snks));
    for n=snks
        %Get response variable and the number of predictors in the stepwise
        %fit. 
        y = autoGLMs{nn}.Variables(:,end);
        nPrdctrs = autoGLMs{nn}.NumPredictors; 
        
        %Build the population that you're sampling from. 
        popToSample = 2:size(X,2);
        respCellInds = find(~cellfun('isempty',regexp(X.Properties.VariableNames,['n',num2str(n),'lag'])));
        popToSample(ismember(popToSample,respCellInds)) = [];
        
        parfor i=1:B
            srcs = [1 randsample(popToSample,nPrdctrs)];

            fitTbl = X(:,srcs);                 %Get design matrix.  
            %fitTbl = X;
            %fitTbl(:,['n',num2str(n)]) = [];    %Get rid of identity.        
            fitTbl(:,end+1) = y;         %Add in response variable. 

            %Do the fit. 
            [~,randcells{i,n}] = SourceSinkGLM(fitTbl,tracetype,false);

            
        end
        
        p.progress;
    end
    p.stop;
end