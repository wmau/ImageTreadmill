TT = []; II = []; IND = [];
for i=243:246
    [t,ind] = getTimePeak(MD(i)); 
    cd(MD(i).Location);
    load(fullfile(pwd,'TimeCells.mat'),'TimeCells');
    load(fullfile(pwd,'TemporalInfo.mat'),'I');
    TT = [TT; t(TimeCells)];
    II = [II; I(TimeCells)];
    IND = [IND; ind(TimeCells)];
end

TTx = [0.25:0.25:10];
m = accumarray(IND,II,[length(TTx),1],@mean,nan);
s = accumarray(IND,II,[length(TTx),1],@std,nan); 
n = sqrt(histc(TT,TTx));
s = s./n;
figure; hold on;
scatter(TT,II,80,'.'); 
errorbar(TTx,m,s,'rd-','linewidth',2)
xlabel('Peak Time [s]'); 
ylabel('Temporal Information [bits/s]');
axis tight; 
set(gca,'ticklength',[0 0]);

[r,p] = corr(TT,II)