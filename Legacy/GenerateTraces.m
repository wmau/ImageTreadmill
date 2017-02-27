load('FinalTraces.mat','rawtrace');

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

minproj = imread('ICmovie_min_proj.tif');
figure;
imshow(minproj,[]);
hold on;
for i=1:n
    PlotNeurons(MD(215),TimeCells(i),c(i,:),3);
end