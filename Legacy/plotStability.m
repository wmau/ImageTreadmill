place = multiProportionStable(MD(215),MD(215),MD(215:2:221),'place');
time = multiProportionStable(MD(215),MD(215),MD(215:2:221),'time'); 
sametuning = TimeCellRemapRate(MD(215),MD(215),MD(215:2:221),[10 10 10 10 10]);

pStable = zeros(1,4);
for i=1:4
    pStable(i) = sum(sametuning(:,i+1)==1)/sum(~isnan(sametuning(:,i+1)));
end

figure;
hold on;
plot(datenum({MD(215:2:221).Date}),place,'-ob','linewidth',3); 
plot(datenum({MD(215:2:221).Date}),time,'-or','linewidth',3); 
legend('Place Cells','Time Cells'); 
datetick('x','mm-dd-yyyy');
ylim([0,1]);
set(gca,'ticklength',[0 0]); 
title('Correlation');
xlabel('Dates'); 
ylabel('Proportion Stable');

figure;
plot(datenum({MD(215:2:221).Date}),pStable,'-ob','linewidth',3); 
datetick('x','mm-dd-yyyy');
ylim([0,1]);
set(gca,'ticklength',[0 0]); 
title('Center of Mass Shift Metric'); 
xlabel('Dates'); 
ylabel('Proportion Stable');