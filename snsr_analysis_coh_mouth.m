%% Loading files
cts_config
setup_environment
group_coh_file = '/media/raphaeltrouille/Sauvegarde/PhD/derivatives/group_cohgrad_all_subjects.mat';
group_coh_perm_file = '/media/raphaeltrouille/Sauvegarde/PhD/derivatives/group_cohgrad_all_subjects_perm.mat';
data_behav = '/mnt/usb-HardDrive_MB/expe_SpeechTrack/analysis/data_behav_dev.mat';
output_folder = '/media/raphaeltrouille/Sauvegarde/PhD/derivatives';
meg_dir = '/mnt/usb-HardDrive_MB/expe_SpeechTrack/meg';

load(group_coh_file)
load(group_coh_perm_file)
load(data_behav)

%% Sensors selection
meg_file = '/mnt/usb-HardDrive_MB/expe_SpeechTrack/meg/meg_2240/meg_2240_set1_order1_vid9.fif';
raw = fiff_setup_read_raw(meg_file, 1);
sensors = get_sensors(raw);
Locc_sensors = sensors.Locc.grads;
Rocc_sensors = sensors.Rocc.grads;

[~, Locc_pair_idx1] = ismember(Locc_sensors, g1);
[~, Locc_pair_idx2] = ismember(Locc_sensors, g2);
Locc_idx = max(Locc_pair_idx1, Locc_pair_idx2);
Locc_paired_sensors = unique(Locc_idx(Locc_idx > 0));


[~, Rocc_pair_idx1] = ismember(Rocc_sensors, g1);
[~, Rocc_pair_idx2] = ismember(Rocc_sensors, g2);
Rocc_idx = max(Rocc_pair_idx1, Rocc_pair_idx2);
Rocc_paired_sensors = unique(Rocc_idx(Rocc_idx > 0));

occipital_pair_indices = [Locc_paired_sensors(:); Rocc_paired_sensors(:)];

%% Params
ref = 4;
f_band = [2, 8]; % Hz
vids_conds = [2 4 6 8]; 
F = find(freq_axis >= f_band(1), 1, 'first') : find(freq_axis <= f_band(2), 1, 'last');
age = data_behav(:, 1);

%% Computation ROI occipital R&L
coh_occ = squeeze(mean(all_cohgrad(:, occipital_pair_indices, F, vids_conds, ref), [3 4]));

%% 1. Global ROI Correlation (Averaging all occipital channels together)
% Average across the sensor dimension (dim 2)
coh_occ_global = mean(coh_occ, 2); 

% Ensure both are column vectors [Subjects x 1]
coh_occ_global = coh_occ_global(:);
age = age(:);

[R_global, p_global] = corr(coh_occ_global, age, 'Type', 'Pearson', 'Rows', 'complete');

fprintf('--- GLOBAL ROI CORRELATION RESULTS ---\n');
fprintf('Global R = %.3f | p-value = %.4f\n\n', R_global, p_global);

%% 2. Sensor-by-Sensor Correlation (Mapping spatial specifics)
[R_sensors, p_sensors] = corr(coh_occ, age, 'Type', 'Pearson', 'Rows', 'complete');

fprintf('--- SIGNIFICANT SENSORS CORRELATION RESULTS ---\n');
for idx = 1:length(occipital_pair_indices)
    if p_sensors(idx) < 0.05
        fprintf('Sensor Pair #%d (Global Matrix Index: %d) : R = %.3f, p = %.4f (*Significant*)\n', ...
            idx, occipital_pair_indices(idx), R_sensors(idx), p_sensors(idx));
    end
end

%% 3. Plotting Global ROI Correlation
figure;
scatter(age, coh_occ_global, 60, 'filled', 'MarkerFaceColor', [0.1216, 0.4706, 0.7059]);
hold on;
lsline; % Add trendline
grid on;
xlabel('Subject Age');
ylabel('Mean Delta/Theta Coherence (Occipital ROI, Video)');
title(sprintf('Coherence vs Age Correlation (R = %.3f, p = %.4f)', R_global, p_global));

%% Topoplot
a = age;

age_ranges = [0 7 8.5 11.5 18 50];
age_ranges(1) = floor(min(a));
age_ranges(end) = ceil(max(a));

group_names = {...
    sprintf('%d-7', age_ranges(1)), ...
    '7-8.5', ...
    '8.5-11.5', ...
    '11.5-18', ...
    '>18'};

n_groups  = length(age_ranges)-1;
n_sensors = size(all_cohgrad,2);

mean_group = nan(n_groups,n_sensors);

for g = 1:n_groups
    if g < n_groups
        idx = a >= age_ranges(g) & a < age_ranges(g+1);
    else
        idx = a >= age_ranges(g) & a <= age_ranges(g+1);
    end
    
    % 1. On extrait la sous-portion (sujets du groupe, tous les capteurs, fréquences F, conditions, ref 4)
    % Taille obtenue : [N_sujets, n_sensors, length(F), 4]
    sub_data = all_cohgrad(idx, :, F, [1 3 5 7], 4); 
    
    % 2. On moyenne sur les dimensions 1 (sujets), 3 (fréquences F) et 4 (conditions)
    % Le résultat aura pour taille : [1, n_sensors, 1, 1]
    mean_sensor_data = mean(sub_data, [1 3 4], 'omitnan');
    
    % 3. On force le passage en vecteur ligne [1 x n_sensors] pour l'assigner proprement
    mean_group(g, :) = mean_sensor_data(:)'; 
    
    fprintf('%s : %d subjects\n', ...
        group_names{g}, sum(idx));
end
size(mean_group)

name1 = cm_labels(g1);
name2 = cm_labels(g2);
cmb_labels = cell(102,1);
for k = 1:102
    n1 = name1{k};
    n2 = name2{k};
    cmb_labels{k} = sprintf('%s+%s', n1, n2(end-3:end));
end

% --- 2. Configuration FieldTrip commune ---
cfg_ft = [];
cfg_ft.xparam   = 'time';
cfg_ft.marker   = 'on';
cfg_ft.style    = 'both';
cfg_ft.colorbar = 'yes';
cfg_ft.layout   = 'neuromag306cmb.lay';
cfg_ft.comment     = 'no';
cfg_ft.colorbar    = 'no';

% (Optionnel) Si vous avez des capteurs significatifs par groupe, 
% il faudra adapter cette section pour chaque groupe.
% cfg_ft.highlight        = 'on';
% cfg_ft.highlightchannel = cmb_labels(significant_sensors); 
% cfg_ft.highlightsymbol  = '*';
% ...

% --- 3. Initialisation de la structure FieldTrip "template" ---
ave_stat = [];
ave_stat.fsample = 1000;
ave_stat.time    = 0;
ave_stat.dimord  = 'chan_time';
ave_stat.label   = cmb_labels;

% --- 4. Boucle de tracé (Option subplot pour tout avoir sur une figure) ---
figure('Position', [100, 100, 1500, 400]); % Figure large pour aligner les groupes

% Trouver la valeur max globale pour avoir la même échelle de couleur (clim) sur tous les plots
max_val = max(mean_group(:));
min_val = min(mean_group(:));

for g = 1:n_groups
    % On isole les données du groupe g (on transpose pour avoir du [102 x 1])
    ave_stat.avg = mean_group(g, :)'; 
    
    % On crée un sous-graphique (1 ligne, 5 colonnes)
    subplot(1, n_groups, g);
    
    % Tracé FieldTrip
    ft_topoplotER(cfg_ft, ave_stat);
    
    % Ajustement pour chaque subplot
    title(group_names{g}, 'FontSize', 12);
    clim([min_val max_val]); % Échelle de couleur uniforme
    colormap(hot);
    hp = get(subplot(1, n_groups, n_groups), 'Position');
    colorbar('Position', [hp(1)+hp(3)+0.01, hp(2), 0.02, hp(4)]);
end

%%
meg_file = '/mnt/usb-HardDrive_MB/expe_SpeechTrack/meg/meg_2240/meg_2240_set1_order1_vid9.fif';
raw = fiff_setup_read_raw(meg_file, 1);
sensors = get_sensors(raw);
Locc_sensors = sensors.Locc.grads;
Rocc_sensors = sensors.Rocc.grads;

[~, Locc_pair_idx1] = ismember(Locc_sensors, g1);
[~, Locc_pair_idx2] = ismember(Locc_sensors, g2);
Locc_idx = max(Locc_pair_idx1, Locc_pair_idx2);
Locc_paired_sensors = unique(Locc_idx(Locc_idx > 0));


[~, Rocc_pair_idx1] = ismember(Rocc_sensors, g1);
[~, Rocc_pair_idx2] = ismember(Rocc_sensors, g2);
Rocc_idx = max(Rocc_pair_idx1, Rocc_pair_idx2);
Rocc_paired_sensors = unique(Rocc_idx(Rocc_idx > 0));

occipital_pair_indices = [Locc_paired_sensors(:); Rocc_paired_sensors(:)];

%% Params
ref = 4;
f_band = [2, 5]; % Hz
vids_conds = [2 4 6 8]; 
F = find(freq_axis >= f_band(1), 1, 'first') : find(freq_axis <= f_band(2), 1, 'last');
age = data_behav(:, 1);

%% Computation ROI occipital R&L
coh_occ = squeeze(mean(all_cohgrad(:, occipital_pair_indices, F, vids_conds, ref), [3 4]));

%% 1. Global ROI Correlation (Averaging all occipital channels together)
% Average across the sensor dimension (dim 2)
coh_occ_global = mean(coh_occ, 2); 

% Ensure both are column vectors [Subjects x 1]
coh_occ_global = coh_occ_global(:);
age = age(:);

[R_global, p_global] = corr(coh_occ_global, age, 'Type', 'Pearson', 'Rows', 'complete');

fprintf('--- GLOBAL ROI CORRELATION RESULTS ---\n');
fprintf('Global R = %.3f | p-value = %.4f\n\n', R_global, p_global);

%% 2. Sensor-by-Sensor Correlation (Mapping spatial specifics)
[R_sensors, p_sensors] = corr(coh_occ, age, 'Type', 'Pearson', 'Rows', 'complete');

fprintf('--- SIGNIFICANT SENSORS CORRELATION RESULTS ---\n');
keep = [];
for idx = 1:length(occipital_pair_indices)
    if p_sensors(idx) < 0.05
        fprintf('Sensor Pair #%d (Global Matrix Index: %d) : R = %.3f, p = %.4f (*Significant*)\n', ...
            idx, occipital_pair_indices(idx), R_sensors(idx), p_sensors(idx));
        keep = cat(2, keep, idx);
    end
end

% --- Plot coh = f(age) for significant sensors
% --- 1. Définir les groupes d'âge ---
% Exemple : Groupes de 10 ans de 20 à 60 ans [20-30[, [30-40[, etc.
% Ajuste ces valeurs selon la répartition réelle de tes sujets !
age_edges = [0 7 8.5 11.5 18 50]; 
num_groups = length(age_edges) - 1;

% Initialisation des matrices pour stocker les résultats
mean_coh = zeros(num_groups, length(keep));
std_coh  = zeros(num_groups, length(keep));
group_labels = cell(1, num_groups);

% --- 2. Calculer la moyenne et l'écart-type par groupe ---
for g = 1:num_groups
    % Trouver les indices des sujets appartenant au groupe g
    idx_subjects = (age >= age_edges(g)) & (age < age_edges(g+1));
    
    if any(idx_subjects)
        % Moyenne et écart-type pour chaque capteur significatif dans ce groupe
        mean_coh(g, :) = mean(coh_occ(idx_subjects, keep), 1, 'omitnan');
        std_coh(g, :)  = std(coh_occ(idx_subjects, keep), 0, 1, 'omitnan');
    else
        mean_coh(g, :) = NaN;
        std_coh(g, :)  = NaN;
    end
    
    % Créer les étiquettes pour l'axe X (ex: "20-30 ans")
    group_labels{g} = sprintf('%d-%d ans', age_edges(g), age_edges(g+1));
end

% --- 3. Représentation Graphique ---
figure;

% Option A : Le diagramme en barres groupées (Idéal si tu as 2 à 4 capteurs max)
b = bar(mean_coh, 'grouped');
hold on;

% Calculer et afficher les barres d'erreur sur chaque barre correspondante
% (Requis pour la rigueur statistique de ton graphique !)
num_bars = size(mean_coh, 2);
for i = 1:num_bars
    % Récupérer les coordonnées X réelles des barres générées par MATLAB
    x_coords = b(i).XEndPoints;
    errorbar(x_coords, mean_coh(:, i), std_coh(:, i), 'k.', 'LineWidth', 1.2);
end

% Esthétique du graphique
set(gca, 'XTickLabel', group_labels);
xlabel('Groupes d''âge');
ylabel('Cohérence Moyenne (± SD)');
title('Cohérence moyenne par groupe d''âge pour les capteurs significatifs');
grid on;

% Légende avec les vrais indices des capteurs
legend(arrayfun(@(x) sprintf('Sensor %d', x), occipital_pair_indices(keep), 'UniformOutput', false), 'Location', 'best');
hold off;

%% 3. Plotting Global ROI Correlation
figure;
scatter(age, coh_occ_global, 60, 'filled', 'MarkerFaceColor', [0.1216, 0.4706, 0.7059]);
hold on;
lsline; % Add trendline
grid on;
xlabel('Subject Age');
ylabel('Mean Delta/Theta Coherence (Occipital ROI, Video)');
title(sprintf('Coherence vs Age Correlation (R = %.3f, p = %.4f)', R_global, p_global));

%% Topoplot Contrast
a = age;

age_ranges = [0 7 8.5 11.5 18 50];
age_ranges(1) = floor(min(a));
age_ranges(end) = ceil(max(a));

group_names = {...
    sprintf('%d-7', age_ranges(1)), ...
    '7-8.5', ...
    '8.5-11.5', ...
    '11.5-18', ...
    '>18'};

n_groups  = length(age_ranges)-1;
n_sensors = size(all_cohgrad,2);

mean_group_vid = nan(n_groups,n_sensors);
mean_group_novid = nan(n_groups,n_sensors);
for g = 1:n_groups
    if g < n_groups
        idx = a >= age_ranges(g) & a < age_ranges(g+1);
    else
        idx = a >= age_ranges(g) & a <= age_ranges(g+1);
    end
    
    % 1. On extrait la sous-portion (sujets du groupe, tous les capteurs, fréquences F, conditions, ref 4)
    % Taille obtenue : [N_sujets, n_sensors, length(F), 4]
    sub_data_vid = all_cohgrad(idx, :, F, [2 4 6 8], 4); 
    sub_data_novid = all_cohgrad(idx, :, F, [1 3 5 7], 4);
    % 2. On moyenne sur les dimensions 1 (sujets), 3 (fréquences F) et 4 (conditions)
    % Le résultat aura pour taille : [1, n_sensors, 1, 1]
    mean_sensor_vid = mean(sub_data_vid, [1 3 4], 'omitnan');
    mean_sensor_novid = mean(sub_data_novid, [1 3 4], 'omitnan');
    % 3. On force le passage en vecteur ligne [1 x n_sensors] pour l'assigner proprement
    mean_group_vid(g, :) = mean_sensor_vid(:)'; 
    mean_group_novid(g, :) = mean_sensor_novid(:)';
    fprintf('%s : %d subjects\n', ...
        group_names{g}, sum(idx));
end
size(mean_group)

name1 = cm_labels(g1);
name2 = cm_labels(g2);
cmb_labels = cell(102,1);
for k = 1:102
    n1 = name1{k};
    n2 = name2{k};
    cmb_labels{k} = sprintf('%s+%s', n1, n2(end-3:end));
end
contrast = (mean_group_vid - mean_group_novid) ./(mean_group_vid + mean_group_novid);
% --- 2. Configuration FieldTrip commune ---
cfg_ft = [];
cfg_ft.xparam   = 'time';
cfg_ft.marker   = 'on';
cfg_ft.style    = 'both';
cfg_ft.colorbar = 'yes';
cfg_ft.layout   = 'neuromag306cmb.lay';
cfg_ft.comment     = 'no';
cfg_ft.colorbar    = 'no';

% (Optionnel) Si vous avez des capteurs significatifs par groupe, 
% il faudra adapter cette section pour chaque groupe.
% cfg_ft.highlight        = 'on';
% cfg_ft.highlightchannel = cmb_labels(significant_sensors); 
% cfg_ft.highlightsymbol  = '*';
% ...

% --- 3. Initialisation de la structure FieldTrip "template" ---
ave_stat = [];
ave_stat.fsample = 1000;
ave_stat.time    = 0;
ave_stat.dimord  = 'chan_time';
ave_stat.label   = cmb_labels;

% --- 4. Boucle de tracé (Option subplot pour tout avoir sur une figure) ---
figure('Position', [100, 100, 1500, 400]); % Figure large pour aligner les groupes

% Trouver la valeur max globale pour avoir la même échelle de couleur (clim) sur tous les plots
max_val = max(contrast(:));
min_val = min(contrast(:));

for g = 1:n_groups
    % On isole les données du groupe g (on transpose pour avoir du [102 x 1])
    ave_stat.avg = contrast(g, :)'; 
    
    % On crée un sous-graphique (1 ligne, 5 colonnes)
    subplot(1, n_groups, g);
    
    % Tracé FieldTrip
    ft_topoplotER(cfg_ft, ave_stat);
    
    % Ajustement cosmétique pour chaque subplot
    title(group_names{g}, 'FontSize', 12);
    clim([min_val max_val]); % Échelle de couleur uniforme
    colormap(hot);
    hp = get(subplot(1, n_groups, n_groups), 'Position');
    colorbar('Position', [hp(1)+hp(3)+0.01, hp(2), 0.02, hp(4)]);
end

