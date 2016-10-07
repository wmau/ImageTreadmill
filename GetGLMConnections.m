function [A,el] = GetGLMConnections(md,glms,neurons)
%
%
%

%%
    cd(md.Location); 
    load('FinalOutput.mat','FT'); 
    nNeurons = size(FT,1); 
    nGLMs = size(glms,2); 
    expMatch = 'n\d+';
    
    A = zeros(nNeurons);
    for i=1:nGLMs
        [s,e] = regexp(glms{i}.CoefficientNames,expMatch);
        
        nParams = length(s); 
        for p=1:nParams
            if ~isempty(s{p})
                trig = str2double(glms{i}.CoefficientNames{p}((s{p}+1):e{p}));
                targ = neurons(i);

                A(trig,targ) = 1; 
            end
        end
    end
    
    el = adj2edgeL(A); 
    
end