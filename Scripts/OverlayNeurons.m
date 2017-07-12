%Warning! Works only if you specify the cell number in the first session. 

clear
close all;
loadMD;

%Sessions and neuron(s). 
mds = MD(292:295);
neurons = [152 58 357]; 
nNeurons = length(neurons);

mapMD = getMapMD(mds);
map = msMatchCells(mapMD,mds,neurons,false); 

%Colors. 
nSessions = length(mds); 
colors = copper(nSessions);

PlotNeurons(mds(1),neurons,colors(1,:),2); hold on;
for s=2:nSessions
    cd(mds(s).Location);
    load('FinalOutput.mat','NeuronImage');
    
    cd(mds(1).Location);
    load(['RegistrationInfo-',mds(s).Animal,'-',mds(s).Date,'-session',...
        num2str(mds(s).Session),'.mat']);
    
    for n=1:nNeurons
        transNeuron = imwarp_quick(NeuronImage{map(n,s)},RegistrationInfoX); 
        b = bwboundaries(transNeuron,'noholes');
        
        x = b{1}(:,1); 
        y = b{1}(:,2); 
        
        plot(y,x,'color',colors(s,:),'linewidth',2);
    end
end
xStart = min(y); 
yStart = min(x); 
line([xStart xStart+10*1.1],[yStart-10 yStart-10],'linewidth',1,'color','k');