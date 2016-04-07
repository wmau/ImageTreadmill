function [I,tCorr,pCorr] = plotStableRhoTI(mapMD,MD1,MDs)
%
%
%
% 
%%
    nMDs = length(MDs); 
    tCorr = cell(nMDs,1);
    pCorr = cell(nMDs,1);
    sf = 10;
    for i=1:nMDs
        [I,~,~,tCorr{i},pCorr{i}] = plotStableTempInfo(mapMD,MD1,MDs(i),'time',0);
    end
    
    close all; 
    Isize = (sf*I/pi).^2;
    
    stable = false(length(I),nMDs);
    unstable = false(length(I),nMDs);
    figure; hold on;
    for i=1:nMDs
        stable(:,i) = tCorr{i}(:,2) < 0.05 & tCorr{i}(:,1) > 0; 
        unstable(:,i) = tCorr{i}(:,2) > 0.05; 
        nStable = sum(stable(:,i)); 
        nUnstable = sum(unstable(:,i)); 
        
        sJitter = i - (0.1 - 0.2*rand(nStable,1));
        uJitter = i - (0.1 - 0.2*rand(nUnstable,1));
        
        scatter(sJitter,tCorr{i}(stable(:,i),1),Isize(stable(:,i)),'go');
        scatter(uJitter,tCorr{i}(unstable(:,i),1),Isize(unstable(:,i)),'ro');
    end
    hold off;
    ylabel('Correlation Rho'); 
    
    figure('position',[130 230 1400 420]);
    for i=1:nMDs
        [~,sP] = corr(I(stable(:,i)),tCorr{i}(stable(:,i),1));
        [~,uP] = corr(I(unstable(:,i)),tCorr{i}(unstable(:,i),1)); 
        
        subplot(1,nMDs,i); hold on;
        scatter(I(stable(:,i)),tCorr{i}(stable(:,i),1),[],'go');
        scatter(I(unstable(:,i)),tCorr{i}(unstable(:,i),1),[],'ro');

        if i==1, ylabel('Intersession Rho'); end
        xlabel('Temporal Information'); 
        
        legend(['p=',num2str(round(sP,3))],['p=',num2str(round(uP,3))],'location','southeast');
    end
end