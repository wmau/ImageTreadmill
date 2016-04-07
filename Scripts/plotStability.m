place = multiProportionStable(MD(244),MD(243),MD(243:247),'place');
time = multiProportionStable(MD(244),MD(243),MD(243:247),'time'); 

figure;
hold on;
plot(datenum({MD(243:247).Date}),place,'-ob','linewidth',3); 
plot(datenum({MD(243:247).Date}),time,'-or','linewidth',3); 
legend('Place Cells','Time Cells'); 
datetick('x','mm-dd-yyyy');
ylim([0,1]);
set(gca,'ticklength',[0 0]); 
xlabel('Dates'); 
ylabel('Proportion Stable');