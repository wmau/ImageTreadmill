function proportionStable(mapMD,MD1,MD2,neurontype)
%
%
%

%%
    dates = {MD1.Date MD2.Date};
    dateTitles = dates; 
    for s=1:2
        dateTitles{s}(3:3:6) = '-';
    end
    
    neurontype = lower(neurontype); 
    switch neurontype
        case 'time'
            load(fullfile(MD1.Location,'TimeCells.mat'),'TimeCells'); 
            n = length(TimeCells);
            [~,r] = PlaceTimeCorr(mapMD,MD1,MD2,TimeCells); 
        case 'place'
            load(fullfile(MD1.Location,'PlaceMaps.mat'),'pval'); 
            PlaceCells = find(pval > 0.95);
            n = length(PlaceCells); 
            [r,~] = PlaceTimeCorr(mapMD,MD1,MD2,PlaceCells); 
    end
    
    p = sum(r(:,2)<0.05)/n;
    figure;
    pc = pie([p 1-p]); 
        ptext = findobj(pc,'Type','text');
        percent = get(ptext,'String'); 
        label = strcat({'Stable: ','Unstable: '},percent');
        pc(1).FaceColor = [0.7 0.7 0.7];    ptext(1).String = label(1);
        pc(3).FaceColor = 'w';              ptext(2).String = label(2); 
        title([dateTitles{1}, ' vs ', dateTitles{2}]);
end