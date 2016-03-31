place = multiProportionStable(MD(241),MD(240),MD(240:244),'place');
time = multiProportionStable(MD(241),MD(240),MD(240:244),'time'); 

figure;
hold on;
plot(datenum({MD(240:244).Date}),place,'-ob','linewidth',3); 
plot(datenum({MD(240:244).Date}),time,'-or','linewidth',3); 
legend('Place Cells','Time Cells'); 
datetick('x','mm-dd-yyyy');
ylim([0,1]);
set(gca,'ticklength',[0 0]); 
xlabel('Dates'); 
ylabel('Proportion Stable');