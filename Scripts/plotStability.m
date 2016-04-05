place = multiProportionStable(MD(243),MD(242),MD(242:246),'place');
time = multiProportionStable(MD(243),MD(242),MD(242:246),'time'); 

figure;
hold on;
plot(datenum({MD(242:246).Date}),place,'-ob','linewidth',3); 
plot(datenum({MD(242:246).Date}),time,'-or','linewidth',3); 
legend('Place Cells','Time Cells'); 
datetick('x','mm-dd-yyyy');
ylim([0,1]);
set(gca,'ticklength',[0 0]); 
xlabel('Dates'); 
ylabel('Proportion Stable');