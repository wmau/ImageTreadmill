loadMD;
load(fullfile(MD(292).Location,'MapsAndPeaks_PCs.mat'));

bin_type = 'u';

%% take all peak times
peaks_tc = [];
plot_on = 1;
for ii = 1:size(PC_DATA,2)
        peaks_tc = [peaks_tc; PC_DATA(ii).peaks(PC_DATA(ii).PlaceCells)];
end


%% plot hist with 50 bins
figure
t1 = [log10(1):0.001:log10(80)];
i = 50;
[n,pos] = hist(peaks_tc,i);
plot(log10(pos),log10(n),'linewidth',2)
xlabel('peak time [s]')
ylabel('number of cells per bin')
%grid on
% set(gca,'xlim',[-1.1 1.1])
% set(gca,'xtick',[-1 0 1])
% set(gca,'xticklabel',[0.1 1 10])
% set(gca,'ytick',[0 1 2 3])
% set(gca,'yticklabel',[0 10 100 1000])
%hold on, plot(log10(t),log10(t.^-1),'b')
pos(find(n==0))=[]; %remove bins with zero elements
n(find(n==0))=[]; %remove bins with zero elements
mdl = fitlm(log10(pos'),log10(n'),'y ~ x1'); %mdl = fitlm(uut,sst,'y ~ x1-1'); without the constat term
x1 = mdl.Coefficients.Estimate;
SE = mdl.Coefficients.SE;
hold on,plot(t1,t1*x1(2)+x1(1),'r','linewidth',2)
%title(sprintf('%i bins, slope %2.2f, SE %2.2f inter %2.2f',i,x1(2),SE(2),x1(1)))
set(gcf,'color','w')
%export_fig(sprintf('TC_PL_fit_allsessions_hist_log10_log10_%s_50bins_linlabel.pdf',bin_type));
%saveas(gcf,sprintf('TC_PL_fit_allsessions_hist_log10_log10_%s_50bins_linlabel.fig',bin_type));