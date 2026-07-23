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
% On fusionne les résultats (si le capteur est dans g1, son index est > 0, sinon il est dans g2)
Locc_idx = max(Locc_pair_idx1, Locc_pair_idx2);
% On élimine les doublons de paires (car g1 et g2 d'une même paire pointent vers le même index)
% et on enlève les 0 au cas où un capteur n'aurait pas été trouvé
Locc_paired_sensors = unique(Locc_idx(Locc_idx > 0));

% On fait exactement la même chose pour le côté droit (Rocc)
[~, Rocc_pair_idx1] = ismember(Rocc_sensors, g1);
[~, Rocc_pair_idx2] = ismember(Rocc_sensors, g2);
Rocc_idx = max(Rocc_pair_idx1, Rocc_pair_idx2);
Rocc_paired_sensors = unique(Rocc_idx(Rocc_idx > 0));

occipital_pair_indices = [Locc_paired_sensors(:); Rocc_paired_sensors(:)];
%% Phase 0.bis
ref = 4;
f_band = [2, 5]; % Hz
F = find(freq_axis >= f_band(1), 1, 'first') : find(freq_axis <= f_band(2), 1, 'last');

% 1. Calcul de la moyenne sur les fréquences (3) et conditions visuelles (4)
coh_visual = squeeze(mean(all_cohgrad(:, :, F, [1 3 5 7], ref), [3 4]));
coh_avg_sub = squeeze(mean(coh_visual, 1));

% Labels
% Construct labels
name1 = cm_labels(g1);
name2 = cm_labels(g2);
cmb_labels = cell(102,1);
for k = 1:102
    n1 = name1{k};
    n2 = name2{k};
    cmb_labels{k} = sprintf('%s+%s', n1, n2(end-3:end));
end

% 2. Préparation de la structure FieldTrip
ave = [];
ave.fsample = 1000;
ave.time    = 0;
ave.dimord  = 'chan_time';
ave.avg     = coh_avg_sub(:); % <-- CORRIGÉ : utilise ton vecteur calculé
ave.label   = cmb_labels;     % Assure-toi que cmb_labels contient tes paires combinées (grad1 + grad2)

% 3. Structure de configuration : cfg
cfg_ft = [];
cfg_ft.xparam   = 'time';
cfg_ft.marker   = 'on';
cfg_ft.style    = 'both';
cfg_ft.colorbar = 'yes';
cfg_ft.layout   = 'neuromag306cmb.lay'; % Layout combiné pour les graduomètres Elekta/Neuromag

% 4. Affichage du Topoplot
fig = figure('Position', [100, 100, 800, 600]);
ft_topoplotER(cfg_ft, ave);
title(sprintf('Group Visual Coherence (Avg Novid conds)(Theta Band: %.1f - %.1f Hz)', f_band(1), f_band(2)));

% 5. Sauvegarde de la figure
% CORRIGÉ : Utilisation de output_folder défini en Phase 0
fig_output_dir = fullfile(output_folder, 'results', 'figures'); 
if ~exist(fig_output_dir, "dir")
    mkdir(fig_output_dir); 
end

% Sauvegarde propre en PNG haute résolution
print(fig, '-dpng', '-r300', fullfile(fig_output_dir, 'group_coherence_theta_novid_conds.png'));
fprintf('\n [OK] Topoplot saved successfully in: %s\n', fig_output_dir);


%% Phase 0
% =========================================================================
% Extract age and subjects IDs
age = data_behav(:, 1);
subjects = table(subject_ids, age);

% -----------
% Correlation subject_nave avec l'age
% Vérifier si nave est une covariable -> moins d'epochs pour les plus
% jeunes ? 

%       per-condition
p_cond = zeros(1, 10);
r_cond = zeros(1, 10);
for cond = 1:10
    [r_cond(cond), p_cond(cond)] = corr(age, subject_nave(:, cond), 'type', 'Spearman');
end
%       Averaged conditions
grand_nave = mean(subject_nave, 2);
[r_grand, p_grand] = corr(age, grand_nave, 'type', 'Spearman');
%  /!\ -> Conclusion :  Pas de significativité (conds/average).

% -----------
% Verification equivalence cond 9 et 10
ref = 4;
f_band = [0.5, 2]; % Hz
F = find(freq_axis >= f_band(1), 1, 'first') : find(freq_axis <= f_band(2), 1, 'last');

coh_9  = squeeze(mean(mean(all_cohgrad(:, :, F,  9, ref), 3), 2));
coh_10 = squeeze(mean(mean(all_cohgrad(:, :, F, 10, ref), 3), 2));

p = signrank(coh_9, coh_10);
% /!\ -> Conclusion: p= 0.7 9 et 10 equivalent au niveau macroscopique

%% Phase 1
% =========================================================================
% REFERENCES| 1: 'global sound'  | 
%           | 2: 'attended sound'|
%           | 3: 'noise sound'   | 
%           | 4: 'mouth opening' |
%
% CONDITIONS | |1|2|3|4|5|6|7|8|9|10 <-
% vid        | [0|1|0|1|0|1|0|1|0| 0;
% energetic  | [0|0|1|1|0|0|1|1|0| 0; <- moyenner 13 56 78
% info       | [0|0|0|0|1|1|1|1|0| 0;
% SiN        | [1|1|1|1|1|1|1|1|0| 0;

% Indice de bénéfice audiovisuel (AVB)
% pairs (1,2) | (3,4) | (5,6) | (7,8)
% Identifier l'apport visuel pour même type de bruit
% -------------------------------------------------------------------------
ref = 4;
f_band = [0.5, 2];
F = find(freq_axis >= f_band(1), 1, 'first') : find(freq_axis <= f_band(2), 1, 'last');
pairs = [1, 2;
         3, 4;
         5, 6;
         7, 8];
normalized_difference = zeros(size(all_cohgrad, 1), length(occipital_pair_indices), size(pairs, 1));
th = mean(all_cohgrad(:) / 10);
% Compute normalized difference (AVB) for each pair Moyenner puis contraste
for i = 1:length(pairs)
    cond_vid   = pairs(i, 2);
    cond_novid = pairs(i, 1);
    coh_vid   = squeeze(mean(all_cohgrad(:, occipital_pair_indices, F, cond_vid, ref), 3));
    coh_novid = squeeze(mean(all_cohgrad(:, occipital_pair_indices, F, cond_novid, ref), 3));
    normalized_difference(:, :, i) = (coh_vid - coh_novid) ./ (coh_vid + coh_novid + th);
end

global_avb = mean(normalized_difference, [2 3]);

[r_contraste, p_contraste] = corr(global_avb, age, 'Type', 'Spearman', 'Rows', 'complete');
for p = 1:length(p_contraste)
    if p_contraste(p) < 0.05
        fprintf("idx=%d/%d -> p=%.3f\n", p, length(p_contraste), p_contraste(p));
    end
end
% -------------------------------------------------------------------------
% Indice de sensibilité au bruit (NSI)
% (silence - bruit) / (silence + bruit)
vid_conds       = [2, 4, 6, 8];
no_vid_conds    = [1, 3, 5, 7];
noiseless_conds = [9, 10];
f_band          = [2, 5];
ref             = 4;
F = find(freq_axis >= f_band(1), 1, 'first') : find(freq_axis <= f_band(2), 1, 'last');


coh_noisy_vid   = squeeze(mean(mean(all_cohgrad(:,occipital_pair_indices, F, vid_conds, ref), 3), 4));
coh_noisy_novid = squeeze(mean(mean(all_cohgrad(:,occipital_pair_indices, F, no_vid_conds, ref), 3), 4));
coh_noiseless   = squeeze(mean(mean(all_cohgrad(:,occipital_pair_indices, F, noiseless_conds, ref), 3), 4));

nsi_vid    = (coh_noiseless - coh_noisy_vid) ./ (coh_noiseless + coh_noisy_vid);
nsi_no_vid = (coh_noiseless - coh_noisy_novid) ./ (coh_noiseless + coh_noisy_novid);



%% 3. Calcul de la corrélation
% On utilise 'Rows','complete' au cas où il y aurait des variables manquantes (NaN)
age = data_behav(:,1);
[R, pV] = corr(coh_noisy_vid, age, 'Type', 'Pearson', 'rows','complete');

fprintf('--- Résultats de la Corrélation ---\n');
fprintf('Coefficient de corrélation (R) : %.3f\n', R);
fprintf('p-value : %.4f\n', pV);

if pV < 0.05
    fprintf('La corrélation est statistiquement significative !\n');
else
    fprintf('Pas de corrélation statistiquement significative.\n');
end

%% 4. Graphique (Scatter Plot)
figure;
scatter(age, coh_noisy_vid, 60, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);
hold on;
% Ajout de la ligne de tendance
lsline; 
grid on;
xlabel('Âge des sujets');
ylabel('Cohérence Delta Moyenne (Occipital - Vidéo)');
title(sprintf('Corrélation entre la Cohérence et l''Âge (R = %.2f, p = %.3f)', R, pV));