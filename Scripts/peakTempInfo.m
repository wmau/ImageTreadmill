TT = []; II = []; IND = [];
for i=240:244
    [t,ind] = getTimePeak(MD(i)); 
    cd(MD(i).Location);
    load(fullfile(pwd,'TimeCells.mat'),'TimeCells');
    load(fullfile(pwd,'TemporalInfo.mat'),'I');
    TT = [TT; t(TimeCells)];
    II = [II; I(TimeCells)];
    IND = [IND; ind(TimeCells)];
end

m = accumarray(IND,II,[],@mean);
s = accumarray(IND,II,[],@std); 
tbl = tabulate(TT);
n = sqrt(tbl(:,2));
s = 2.*s./n;
figure; hold on;
scatter(TT,II,80,'.'); 
errorbar(unique(TT),m,s,'rd-','linewidth',2)
xlabel('Peak Time [s]'); 
ylabel('Temporal Information [bits/s]');
axis tight; 
set(gca,'ticklength',[0 0]);

[r,p] = corr(TT,II)