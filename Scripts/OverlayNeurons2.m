s1 = MD(215);
s2 = MD(217);

load(fullfile(s1.Location,...
    'RegistrationInfo-GCamp6f_45_treadmill-12_01_2015-session4.mat'));

load(fullfile(s1.Location,'FinalOutput.mat'),'FT');
n1 = size(FT,1); 
load(fullfile(s2.Location,'FinalOutput.mat'),'NeuronImage'); 
n2 = size(NeuronImage,2); 

for n=1:n2
    transNeuron = imwarp_quick(NeuronImage{n},RegistrationInfoX);
    b = bwboundaries(transNeuron,'noholes'); 
    
    x{n} = b{1}(:,1);
    y{n} = b{1}(:,2); 
end

PlotNeurons(s1,1:n1,[0 0 1 .4],2); 
hold on;
for n=1:n2
    plot(y{n},x{n},'color',[1 .65 0 .4],'linewidth',2);
end