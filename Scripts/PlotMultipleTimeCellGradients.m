session1 = MD(292);
session2 = MD(293); 

[~,peak1] = getTimePeak(session1); 
dotSize1 = peak1*50;
[~,peak2] = getTimePeak(session2);
dotSize2 = peak2*50;

matchedCells = msMatchMultiSessionCells([session1 session2],...
    {getTimeCells(session1) getTimeCells(session2)}); 
regInfo = image_registerX(session1.Animal,session1.Date,session1.Session,...
    session2.Date,session2.Session);

day1cells = matchedCells(matchedCells(:,1)>0,1);
day2cells = matchedCells(matchedCells(:,2)>0,2);

figure('Position',[460 125 990 810]);
PlotNeuronsDots(session1,'neurons',day1cells,'dotSizes',...
    dotSize1(day1cells),'dotColors','b','filled',true);
figure('Position',[460 125 990 810]);
PlotNeuronsDots(session2,'neurons',day2cells,'dotSizes',...
     dotSize2(day2cells),'dotColors','r','filled',true,...
        'transformation',regInfo); 

figure('Position',[460 125 990 810]);
hold on;
PlotNeuronsDots(session1,'neurons',day1cells,'dotSizes',...
    dotSize1(day1cells),'dotColors','b','filled',true,'transparency',0.4);
PlotNeuronsDots(session2,'neurons',day2cells,'dotSizes',...
    dotSize2(day2cells),'dotColors','r','filled',true,...
     'transparency',0.4,'transformation',regInfo); 
% figure;
% difference = abs((peak1(matchedCells(:,1)) - peak2(matchedCells(:,2))))*20+0.01;
% PlotNeuronsDots(MD(292),'neurons',matchedCells(:,1),'dotSizes',difference);

s1proj = imread(fullfile(session1.Location,'ICmovie_min_proj.tif')); 
s2proj = imread(fullfile(session2.Location,'ICmovie_min_proj.tif')); 
s2proj = imwarp_quick(s2proj,regInfo); 

figure(51);
imagesc(s1proj); colormap gray; axis equal; axis off;
figure(52);
imagesc(s2proj); colormap gray; axis equal; axis off;
figure(50);
imagesc(s1proj-s2proj); colormap gray; axis equal; axis off; 


minGray = min([min(s1proj) min(s2proj)]); 
maxGray = max([max(s1proj) max(s2proj)]); 
figure(51); caxis([minGray maxGray]); 
figure(52); caxis([minGray maxGray]); 