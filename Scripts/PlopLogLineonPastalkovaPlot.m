%% Just for putting a scalar fit line on a Pastalkova plot. 

set(gcf,'position',[680 370 320 610]); 
Nts = size(ensemble,1); %number of cells
min_t = 200; %ms
max_t = 10000; %ms
ts = exp(linspace(log(min_t),log(max_t),Nts)); % peak times
hold on ,plot(ts/1000,1:length(ts),'r--', 'linewidth',3)
set(gca,'fontsize',15)