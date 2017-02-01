clear
close all;
loadMD;
s1 = MD(292);
s2 = MD(293);

cd(s1.Location);
mapMD = getMapMD([s1,s2]);

load(fullfile(s1.Location,...
    'RegistrationInfo-GCamp6f_45_treadmill-12_01_2015-session10.mat'));

load(fullfile(s1.Location,'FinalOutput.mat'),'PSAbool');
n1 = size(PSAbool,1); 
load(fullfile(s2.Location,'FinalOutput.mat'),'NeuronImage'); 
n2 = size(NeuronImage,2); 

for n=1:n2
    transNeuron = imwarp_quick(NeuronImage{n},RegistrationInfoX);
    b = bwboundaries(transNeuron,'noholes'); 
    
    try
        x{n} = b{1}(:,1);
        y{n} = b{1}(:,2); 
    catch
        x{n} = nan;
        y{n} = nan;
    end
end

PlotNeurons(s1,1:n1,[0 0 1 .4],2); 
hold on;
for n=1:n2
    plot(y{n},x{n},'color',[1 0 0 .4],'linewidth',2);
end

line([0 100*1.1],[0 0],'linewidth',5,'color','k');