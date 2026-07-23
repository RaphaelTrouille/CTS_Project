% Loading files
group_coh_file = '/media/raphaeltrouille/Sauvegarde/PhD/derivatives/group_cohgrad_all_subjects.mat';
group_coh_perm_file = '/media/raphaeltrouille/Sauvegarde/PhD/derivatives/group_cohgrad_all_subjects_perm.mat';
data_behav = '/mnt/usb-HardDrive_MB/expe_SpeechTrack/analysis/data_behav_dev.mat';
output_folder = '/media/raphaeltrouille/Sauvegarde/PhD/derivatives';

load(group_coh_file)
load(group_coh_perm_file)
load(data_behav)


% (subjects x sensors x frequencies x conditions x references)

% -------------------------------------------------------------------------
% Frequencies of interest:
f_band = [0.5 2];
Find_start = find(freq_axis >= f_band(1), 1, 'first');
Find_end   = find(freq_axis <= f_band(2), 1, 'last');
F = Find_start:Find_end;

fprintf("Selected band: : [%.2f - %.2f]\n", freq_axis(Find_start), freq_axis(Find_end));

% -------------------------------------------------------------------------
% Sorting by ages:
age_ranges = [0 7 8.5 11.5 18 50];
subjects_age = data_behav(:,1);
[sub_ages_sorted, sort_idx] = sort(subjects_age, 'ascend');
age_sorted_coh = all_cohgrad(sort_idx, :, :, :, :);
fprintf("Youngest: %.2f | Oldest: %.2f | Mean age: %.2f | std: %.2f\n", ...
    sub_ages_sorted(1), sub_ages_sorted(end), mean(sub_ages_sorted), std(sub_ages_sorted));
all_cohgrad = age_sorted_coh;
% -------------------------------------------------------------------------
% Extraction & computation of the noise resistance index
% REFERENCES| 1: 'global sound'  | 
%           | 2: 'attended sound'|
%           | 3: 'noise sound'   | 
%           | 4: 'mouth opening' |
%
% CONDITIONS | |1|2|3|4|5|6|7|8|9|
%            | |_|_|_|_|_|_|_|_|_|
% vid        | [0|1|0|1|0|1|0|1|0];
% energetic  | [0|0|1|1|0|0|1|1|0];
% info       | [0|0|0|0|1|1|1|1|0];
% SiN        | [1|1|1|1|1|1|1|1|0];
n_ref = 4;

% 1. Average over frequency band of interest
coh_f_ref = squeeze(mean(all_cohgrad(:, :, F, :, n_ref), 3));

% 2. Separate and average across conditions
coh_noisy = squeeze(mean(coh_f_ref(:, :, [2 4 6 8]), 3));
coh_noiseless = squeeze(mean(coh_f_ref(:, :, 9:10), 3)); %<- 1 3 5 7

% 3. Apply baseline floor threshold
th = mean(all_cohgrad(:)) / 10;
coh_noisy_th     = max(coh_noisy, th);
coh_noiseless_th = max(coh_noiseless, th);

% 4. Normalized contrast (Noise Resistance Index)
% Close to 0 -> Poorly sensitive to noise
% Close to 1 -> Highly sensitive to noise
coh_NoiseRes = (coh_noiseless_th - coh_noisy_th) ./ (coh_noiseless_th + coh_noisy_th);
fprintf("\n [OK] 'Coh_NoiseRes' index computed successfully! Matrix size: [%dx%d]\n", ...
    size(coh_NoiseRes, 1), size(coh_NoiseRes, 2));

% -------------------------------------------------------------------------
% Screening: Find sensors which NoiseRes significantly change by the age

% Preallocate space
r_age = zeros(1, size(coh_NoiseRes, 2));
p_age = zeros(1, size(coh_NoiseRes, 2));

% Compute correlation subject by subject for each sensor
for n_src = 1:size(coh_NoiseRes, 2)
    [r_age(n_src), p_age(n_src)] = corr(sub_ages_sorted, coh_NoiseRes(:, n_src), 'type', 'Spearman');
end

% Significant sensors
significant_sensors = find(p_age < 0.05);
fprintf('Number of sensors correlatad to age : %d / %d', ...
    length(significant_sensors), size(coh_NoiseRes, 2));

% Most significant sensor
[min_p, best_sensor] = min(p_age);
figure('Position', [100, 100, 800, 500]);

scatter(sub_ages_sorted, coh_NoiseRes(:, best_sensor), 45, 'filled', ...
    'MarkerFaceColor', [0.1 0.5 0.7], 'MarkerEdgeColor', 'w', 'MarkerFaceAlpha', 0.6);
hold on

% Non-Linar local trend line (LOWESS)
smoothed_trend = smooth(sub_ages_sorted, coh_NoiseRes(:, best_sensor), 0.4, 'lowess');
plot(sub_ages_sorted, smoothed_trend, 'r', LineWidth=3);

grid on;
set(gca, 'GridLineStyle', ':', 'LineWidth', 1.2, 'FontSize', 12);
xlabel('Age of Subjects (years)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Noise Resistance Index (Contrast)', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Developmental Trajectory of Visual Noise Resistance\n(Sensor #%d | p = %.4f)', best_sensor, min_p), ...
      'FontSize', 14, 'FontWeight', 'bold');

legend({'Subjects (N=144)', 'Developmental Trend (LOWESS)'}, 'Location', 'Best');

% -------------------------------------------------------------------------
% Topoplot de r ou p.

name1 = cm_labels(g1);
name2 = cm_labels(g2);
cmb_labels = cell(102,1);
for k = 1:102
    n1 = name1{k};
    n2 = name2{k};
    cmb_labels{k} = sprintf('%s+%s', n1, n2(end-3:end));
end

data_to_plot = -log10(p_age(:));

ave_stat = [];
ave_stat.fsample = 1000;
ave_stat.time    = 0;
ave_stat.dimord  = 'chan_time';
ave_stat.label   = cmb_labels;
ave_stat.avg     = data_to_plot;

cfg_ft = [];
cfg_ft.xparam   = 'time';
cfg_ft.marker   = 'on';
cfg_ft.style    = 'both';
cfg_ft.colorbar = 'yes';
cfg_ft.layout   = 'neuromag306cmb.lay';

cfg_ft.highlight        = 'on';
cfg_ft.highlightchannel = cmb_labels(significant_sensors);
cfg_ft.highlightsymbol  = '*';
cfg_ft.highlightcolor   = [1 1 1];
cfg_ft.highlightsize    = 10;

fig = figure('Position', [100, 100, 800, 600]);
ft_topoplotER(cfg_ft, ave_stat);

clim([0 max(ave_stat.avg)]);
colormap(hot);

title('Statistical significance Map (-log(p))');

if ~exist(output_folder, 'dir'), mkdir(output_folder); end
print('-dpng', fullfile(output_folder, 'topoplot_spearman_age_NoiseRes.png'))





