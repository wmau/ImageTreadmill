clear;
loadMD;

fulldataset = [MD(292:303) MD(305:308)];            %Sessions.
nSessions = length(fulldataset); 

overlap = cell(1,nSessions); 
for s=1:nSessions
    [trial_density_map,neurons,order] = trialDensityMap(fulldataset(s));
    
    nTrials = size(trial_density_map,2); 
    
    [cells,~] = find(trial_density_map(:,1:8));
    original_ensemble = unique(cells); 
    original_ensemble_size = length(original_ensemble); 
    
    overlap{s} = [];
    for t=1:5:nTrials
        try
            [active_cells,~] = find(trial_density_map(:,t:t+4)); 
        catch
            [active_cells,~] = find(trial_density_map(:,t:end)); 
        end
        active_cells = unique(active_cells);
        same_cells = intersect(active_cells, original_ensemble); 
        
        overlap{s} = [overlap{s} length(same_cells)/original_ensemble_size]; 
    end
end

longest_session = max(cellfun('length',overlap)); 
overlap_all = nan(nSessions,longest_session); 

for s=1:nSessions
    nTrials = length(overlap{s});
    
    for t=1:nTrials
        overlap_all(s,t) = overlap{s}(t);
    end
end

overlap_all(:,end-2:end) = [];

overlap_mean = nanmean(overlap_all)*100;
overlap_sem = standarderror(overlap_all)*100;
errorbar(overlap_mean,overlap_sem,'capsize',0);
ax = gca; 
ylabel('% time cells overlapping');
xlabel('Trial block'); 
xlim([0,longest_session-2]);
make_plot_pretty(ax);

X = overlap_all(:);
trial_num = repmat(1:longest_session-3,nSessions,1);
trial_num = trial_num(:);

[p,stats,tbl] = anovan(X,{trial_num});