cd(MD(292).Location);
load('FinalOutput.mat','NeuronTraces');
NumNeurons = size(NeuronTraces.RawTrace,1);
TimeCells = getTimeCells(MD(292));
TimeCells(9:10) = TimeCells(13:14); 
rawtrace = NeuronTraces.RawTrace;
t = (1:size(rawtrace,2))./20;

figure('Position',[590 240 340 640]);
n = 10;
colors = distinguishable_colors(n);
c = zeros(NumNeurons,3);
c(TimeCells(1:n),:) = colors; 
for i=1:n
    AX(i) = subplot(10,1,i);
    plot(t,rawtrace(TimeCells(i),:),'color',c(TimeCells(i),:));
    ylim([min(rawtrace(TimeCells(i),:)), max(rawtrace(TimeCells(i),:))]);
    axis off; 
end

set(AX,'YLim',[min([AX.YLim]),max([AX.YLim])],...
    'XLim',[0,max(t)]);
axis on;

proj = imread('MaxProj.tif');
figure('Position',[950 260 760 620]);
imshow(proj,[]);
hold on;
PlotNeurons(MD(292),TimeCells(1:n)',c,3); 
line([0 100*1.10],[5 5],'linewidth',5,'color','w');