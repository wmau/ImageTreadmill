clear;
loadMD; close all;
cd(MD(292).Location);
load('FinalOutput.mat','NeuronTraces');
TimeCells = getTimeCells(MD(292));
rawtrace = NeuronTraces.RawTrace;
t = (1:size(rawtrace,2))./20;

figure;
n = 10;
c = rand(n,3);
for i=1:n
    AX(i) = subplot(10,1,i);
    plot(t,rawtrace(TimeCells(i),:),'color',c(i,:));
    ylim([min(rawtrace(TimeCells(i),:)), max(rawtrace(TimeCells(i),:))]);
    axis off; 
end

set(AX,'YLim',[min([AX.YLim]),max([AX.YLim])],...
    'XLim',[0,max(t)]);
axis on; set(gca,'ticklength',[0 0]);

movie = 'BPDFF.h5';
% Get movie information. 
[Xdim,Ydim,NumFrames] = Get_T_Params('Xdim','Ydim','NumFrames');

% step 1 build up a maximum projection, using every 5th frame
proj = zeros(Xdim,Ydim);
for i = 1:5:NumFrames
    temp = LoadFrames(movie,i);
    proj(temp > proj) = temp(temp > proj);
end

figure;
imshow(proj,[]);
hold on;
for i=1:n
    PlotNeurons(MD(292),TimeCells(i),c(i,:),3);
end