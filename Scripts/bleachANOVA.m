clear;

loadMD;

%Choose session here. 
sessions = [MD(292) MD(296) MD(300) MD(305)]; 

for s=1:length(sessions)
    session = sessions(s);
    
    %Go to session and load traces. 
    cd(session.Location); 
    load('Pos_align.mat','DFDTtrace'); 

    trace = mean(DFDTtrace);
    T = (length(trace)/20)/60;                  %Time in minutes.
    
    %Bin fluorescence into minute-long bins. 
    [m,sem,chunkID] = bleachCheck(trace,round(T)); 

    subplot(2,2,s);
    errorbar(linspace(0,T,length(m)),m,sem,'capsize',0,'linewidth',4); 
    xlim([0 T]);
    ylim([-.005 .005]);
    xlabel('Time (min)','fontsize',15);
    ylabel('Average delta fluorescence (AU)','fontsize',15);
    set(gca,'tickdir','out','fontsize',12);

    [~,~,stats] = anovan(trace,{chunkID}); 
end