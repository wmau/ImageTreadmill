clear
load('Peaks.mat')
load('Maps.mat')

for an_no = 1:size(MAPS,2) % number of animals
    %compute transfer_map which stores which dates in MAPS correspond to which dates in Peaks
    clear transfer_map 
    for ses_peaks = 1:size(Peaks,2)
        for ses_maps = 1:size(MAPS(an_no).session,2)
            if strcmp(MAPS(an_no).session(ses_maps).Animal,Peaks(ses_peaks).Animal) &&...
                strcmp(MAPS(an_no).session(ses_maps).Date,Peaks(ses_peaks).Date)
                transfer_map(ses_maps) = ses_peaks;
            end          
        end
    end
    % find and aggregate peak times of all the cells
    time_cell_no = 0; % counts how many time cells are there for this animal
    for cell_no = 1:size(MAPS(an_no).map,1)
        time_cell = 0; % tracks if the current cell has ever been recorded as time cell
        no_repeats = 0; % counts number of same cell was recorded
        no_repeats_tc = 0; % number of times same cells was recorded and was clasified as a time cell
        for ses_maps = 1:size(MAPS(an_no).session,2)
            %Peaks(transfer_map(ses_maps)).binned(Peaks(transfer_map(ses_maps)).TimeCells(MAPS(an_no).map(tc_no,ses_maps)))
            if MAPS(an_no).map(cell_no,ses_maps) ~= 0 && ~isnan(MAPS(an_no).map(cell_no,ses_maps))
                no_repeats = no_repeats + 1;
                peak_times{an_no,cell_no,no_repeats} =...
                    Peaks(transfer_map(ses_maps)).binned(MAPS(an_no).map(cell_no,ses_maps));
                if sum(cell_no==Peaks(transfer_map(ses_maps)).TimeCells) %is it labeled as time cell 
                    if ~time_cell
                        time_cell_no = time_cell_no + 1; %if the time cell counter was not already increased
                    end
                    time_cell = 1;
                    no_repeats_tc = no_repeats_tc + 1;
                    peak_times_tc{an_no,time_cell_no,no_repeats_tc} =...
                    Peaks(transfer_map(ses_maps)).binned(MAPS(an_no).map(cell_no,ses_maps));
                else
                    time_cell = 0;
                end
            end
        end
    end
end


for an_no = 1:size(peak_times,1) % number of animals
    for cell_no = 1:size(peak_times,2)
        peak_times_tmp = squeeze(cell2mat(peak_times(an_no,cell_no,:))); 
        if ~isempty(peak_times_tmp)
            peak_times_mean{an_no,cell_no} = mean(peak_times_tmp);
            peak_times_se{an_no,cell_no} = std(peak_times_tmp)/sqrt(length(peak_times_tmp));
        end
    end
    for cell_no = 1:size(peak_times_tc,2)
        peak_times_tmp = squeeze(cell2mat(peak_times_tc(an_no,cell_no,:))); 
        if ~isempty(peak_times_tmp)
            peak_times_tc_mean{an_no,cell_no} = mean(peak_times_tmp);
            peak_times_tc_se{an_no,cell_no} = std(peak_times_tmp)/sqrt(length(peak_times_tmp));
        end
    end
end

figure
for an_no = 1:size(peak_times,1) % number of animals     
    [peak_times_mean_sort, peak_times_mean_sort_ind] = sort(cell2mat(peak_times_mean(an_no,:)));
    peak_times_se_tmp = cell2mat(peak_times_se(an_no,:));
    peak_times_se_sort = peak_times_se_tmp(peak_times_mean_sort_ind);
    subplot(size(peak_times,1),1,an_no), errorbar(peak_times_mean_sort,peak_times_se_sort)
    ylabel(sprintf('Animal %i',an_no))
end
set(gcf,'color','w')
xlabel('Cell #')
export_fig(sprintf('results_PL_fit/SE_peak_time_all_cells.pdf'));
saveas(gcf,sprintf('results_PL_fit/SE_peak_time_all_cells.fig'));


figure
for an_no = 1:size(peak_times,1) % number of animals     
    [peak_times_mean_sort, peak_times_mean_sort_ind] = sort(cell2mat(peak_times_tc_mean(an_no,:)));
    peak_times_se_tmp = cell2mat(peak_times_tc_se(an_no,:));
    peak_times_se_sort = peak_times_se_tmp(peak_times_mean_sort_ind);
    subplot(size(peak_times,1),1,an_no), errorbar(peak_times_mean_sort,peak_times_se_sort)
    ylabel(sprintf('Animal %i',an_no))
end
set(gcf,'color','w')
xlabel('Cell #')
export_fig(sprintf('results_PL_fit/SE_peak_time_time_cells.pdf'));
saveas(gcf,sprintf('results_PL_fit/SE_peak_time_time_cells.fig'));

save('peak_times_tc_mean.mat','peak_times_tc_mean')

peaks_min = 0;
peaks_max = 0;
plot_on = 1;
for an_no = 1:size(peak_times,1)
    [ml_exp(an_no,:), ml_confint_low(an_no,:), ml_confint_high(an_no,:)] =...
        fit_power_law_to_peak_times(cell2mat(peak_times_tc_mean(an_no,:)),peaks_min,peaks_max,plot_on);
    export_fig(sprintf('results_PL_fit/PL_fit_animal_%i.pdf',an_no));
    saveas(gcf,sprintf('results_PL_fit/PL_fit_animal_%i.fig',an_no));
end

figure,errorbar(1:4, ml_exp, ml_confint_low-ml_exp, ml_confint_high-ml_exp)
hold on,plot([0 5],[1 1],'r')
set(gca,'xlim', [0 5])
set(gca,'xtick',[1:4])
xlabel('Animal')
ylabel('-Exponent')
set(gcf,'color','w')
title('Averaged peak times of time cells')
export_fig(sprintf('results_PL_fit/TC_nooverlap_PL_fit_per_animal.pdf'));
saveas(gcf,sprintf('results_PL_fit/TC_nooverlap_PL_fit_per_animal.fig'));