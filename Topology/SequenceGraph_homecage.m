function A = SequenceGraph_homecage(animal,date,session,window,alpha)
%A = SequenceGraph_homecage(animal,date,session,window,alpha)
%
%   
%

%% Load the traces. 
    ChangeDirectory(animal,date,session);
    load(fullfile(pwd,'ProcOut.mat'),'FT');
    
%% Setup.
    %Useful variables. 
    [nNeurons,nFrames] = size(FT);
    dt = 0.05;                          %seconds.
    nBins = window/dt;          
    [onset,~] = getFEpochs(FT);         %Get onset indices for calcium events. 
    null = randi([-nBins,nBins],1,100000);    %Null distribution for comparing empirical lags. 
    A = zeros(nNeurons); 
    
%% Construct the graph. 
    %For each neuron...
    p = ProgressBar(nNeurons);
    for n=1:nNeurons
        %Look at the onset plus/minus the lag.
        wBack = onset{n}-nBins; 
        wForward = onset{n}+nBins; 
        
        %Number of calcium events. 
        nEpochs = length(onset{n});
        
        %Preallocate. Reset for every neuron. Each cell entry represents
        %the lag at which another neuron was coincident within dt seconds. 
        corrDists = cell(nNeurons,1); 
        corrDists{n} = nan;
        
        %For each calcium event...
        for e=1:nEpochs
            b = wBack(e);       %Index of beginning of epoch, minus lag.
            f = wForward(e);    %Index of end of epoch, plus lag.
            
            %Neurons that fire within the window of neuron n. 
            coincident = find(cellfun(@any,cellfun(@(c) c>b & c<f, onset,'unif',0)));
            coincident(coincident==n) = [];             %Remove self-coincidence.
            
            %For each neuron that fires within the window...
            for i=1:length(coincident)
                c = coincident(i);
                lags = b+nBins-onset{c};
                corrDists{c} = [corrDists{c}, lags(lags<nBins & lags>-nBins)];      
            end
            
        end
        
        %For each neuron that was coincident...
        for i=1:length(coincident)
            c = coincident(i);
            
            if ~isempty(corrDists{c})
                %Test lag distribution against a uniform null. 
                h = kstest2(corrDists{c},null,'alpha',alpha);
                
                %Get the average lag to find whether it is positive or
                %negative. 
                avgLag = mean(corrDists{c}); 
                
                %Directed graph. 
                if avgLag>0 && h
                        A(n,c) = avgLag; 
                        Au(n,c) = avgLag; Au(c,n) = avgLag;
                elseif avgLag<0 && h
                    A(c,n) = -avgLag; 
                    Au(n,c) = -avgLag; Au(c,n) = -avgLag;
                end
            end
        end
        
    p.progress;
    end
    p.stop;
    
    %Save data.
    save(['A',num2str(alpha),'.mat'], 'A','Au','window');
end