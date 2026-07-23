% Init
cfg = cts_config;
[meg_dir, ~, deriv_dir, ~, ~] = setup_environment;

% Parameters
n_ref = 2;
Find = 2;
n_cond = 9;
Nkeep = 24;
matfiles_dir = '/mnt/usb-HardDrive_MB/expe_SpeechTrack/matfiles';
perm_stat = 0;
suffix = iif(perm_stat,'_main_coh_perm_stat.mat', '_main_coh.mat');
num_channels = 102;

n_subjects = length(cfg.subjects.include);
all_coh = nan(length(cfg.subjects.include), num_channels);

% Extracting sensors config
first_subject = cfg.subjects.include{1};
meg_pattern = fullfile(meg_dir, first_subject, '*_tsss_mc.fif');
fif_files = dir(meg_pattern);
if ~isempty(fif_files)
    meg_file = fullfile(meg_dir, first_subject, fif_files(1).name);
    raw = fiff_setup_read_raw(meg_file);
    sensors = get_sensors(raw);
    all_ch_names = {raw.info.chs.ch_name}';

    Locc_names = all_ch_names(sensors.picksMEG(sensors.Locc.grads));
    Rocc_names = all_ch_names(sensors.picksMEG(sensors.Rocc.grads));

    l_occ_indices = sensors.Locc.grads;
    r_occ_indices = sensors.Rocc.grads;
else
    disp("No  *_tsss_mc.fif files for\n");
end



for i = 1:n_subjects    
    subject = cfg.subjects.include{i};

    fprintf('Subject: %s\n', subject);
    matfile = fullfile(matfiles_dir, subject, [subject suffix]);
    
    % Load matfile
    if ~exist(matfile, "file")
        fprintf("File %s not found\n", matfile)
        continue
    end
    load(matfile, 'CMall')
    all_coh(i, :) = CMall(n_cond).cohgrad(:, Find, n_ref)';

end

mean_coh_group = mean(all_coh, 1, 'omitnan');
% Build gradiometers pairs -> Physical sensors
cm_labels = CMall(n_cond).label;
g1 = double(CMall(n_cond).g1);
g2 = double(CMall(n_cond).g2);



l_occ_names = all_ch_names(l_occ_indices);
r_occ_names = all_ch_names(r_occ_indices);

l_occ_indices_cm = find(ismember(cm_labels, l_occ_names));
r_occ_indices_cm = find(ismember(cm_labels, r_occ_names));

% --- Find pairs ---
l_occ_pairs = find(ismember(g1, l_occ_indices_cm) | ismember(g2, l_occ_indices_cm));
r_occ_pairs = find(ismember(g1, r_occ_indices_cm) | ismember(g2, r_occ_indices_cm));

% --- Left hemisphere ---
coh_Locc = mean_coh_group(l_occ_pairs);
[~, sort_idx_L] = sort(coh_Locc, 'descend');
lkeep = l_occ_pairs(sort_idx_L(1:min(Nkeep, length(sort_idx_L))));
lkeep_names = cm_labels(g1(lkeep));

% --- Right hemisphere ---
coh_Rocc = mean_coh_group(r_occ_pairs);
[~, sort_idx_R] = sort(coh_Rocc, 'descend');
rkeep = r_occ_pairs(sort_idx_R(1:min(Nkeep, length(sort_idx_R))));
rkeep_names = cm_labels(g2(rkeep));

% Construct labels
name1 = cm_labels(g1);
name2 = cm_labels(g2);
cmb_labels = cell(102,1);
for k = 1:102
    n1 = name1{k};
    n2 = name2{k};
    cmb_labels{k} = sprintf('%s+%s', n1, n2(end-3:end));
end


















%% ===================================================================
%  ===================================================================
% TOPOPLOT
ave = [];
ave.fsample = 1000;
ave.time    = 0;
ave.dimord  = 'chan_time';
ave.avg     = mean_coh_group(:);
ave.label = cmb_labels;

% Structure preparation: cfg
cfg_ft = [];
cfg_ft.xparam = 'time';
cfg_ft.marker = 'on';
cfg_ft.style  = 'both';
cfg_ft.colorbar = 'yes';
cfg_ft.layout = 'neuromag306cmb.lay';

% Plot
fig = figure('Position', [100, 100, 800, 600]);
ft_topoplotER(cfg_ft, ave);
title(sprintf('Group Coherence (No-Noise - F = %d)', Find));

% Save the figure
fig_output_dir = fullfile(deriv_dir, 'results', 'figures');
if ~exist(fig_output_dir, "dir"), mkdir(fig_output_dir); end
print('-dpng', fullfile(fig_output_dir, 'group_coherence_2Hz.png'))





































