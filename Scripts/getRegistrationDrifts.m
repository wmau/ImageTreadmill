clear
close all;
loadMD; 

s1 = MD(292);
s2 = MD(293); 

mapMD = getMapMD([s1,s2]);

load(fullfile(s1.Location,...
    'RegistrationInfo-GCamp6f_45_treadmill-12_01_2015-session10.mat'));

load(fullfile(s1.Location,'FinalOutput.mat'),'NumNeurons','NeuronImage');
n1 = NeuronImage;
load(fullfile(s2.Location,'FinalOutput.mat'),'NeuronImage');
n2 = NeuronImage;

matches = msMatchCells(mapMD,MD(292:293),1:NumNeurons,true);
N = size(matches,1);

BinBlobs_reg{1} = n1(matches(:,1));

for n=1:N
    BinBlobs_reg{2}{n} = imwarp_quick(n2{matches(n,2)},RegistrationInfoX);
end

[~,d,~,r,~,o] = dist_bw_reg_sessions(BinBlobs_reg,0);

figure;
histogram(d(:,2,:),40,'normalization','probability','edgecolor','none');
xlabel('Centroid distance [\mum]');
ylabel('Proportion');
set(gca,'tickdir','out');

figure;
histogram(o(:,2,:),40,'normalization','probability','edgecolor','none');
xlabel('Orientation difference [degrees]');
ylabel('Proportion');
set(gca,'tickdir','out');

figure;
A = zeros(size(BinBlobs_reg{1}{1}));
for n=1:N 
    A = A + BinBlobs_reg{1}{n} + BinBlobs_reg{2}{n};
end
A(A>2) = 2; 
imagesc(A);
hold on;
axis equal; axis off; 
line([50 150],[500 500],'linewidth',5,'color','k');