%% s07_coh_age_analysis.m
%
% 2-5Hz source-space coherence with the mouth-opening reference, pooled 
% accross condition sets (1,3,5,7) = "novid" and (2,4,6,8)="vid", observed
% vs. surrogate, data-driven occipital ROI, baseline-corrected age
% correlation.
% 
%Decisions locked in:
%   - Reference channel: mouth-opening (n_ref = 4)
%   - ROI: data_driven peak/cluster within an occipital search mask
%   - group-level map: observed vs. surrogate only (per condition set)
%   - age correlation: baseline-corrected coherence (observed - surrogate)
%
% Can be run section by section...

%% 0. Setup
cfg = cts_config;
[meg_dir, ~, deriv_dir, ~, ~] = setup_environment;
exp_dir = '/mnt/usb-HardDrive_MB/expe_SpeechTrack';
subjects_dir = '/mnt/usb-HardDrive_MB/expe_SpeechTrack/subjects';

n_ref = 4;  % mouth opening
subjects = cfg.subjects.include;    % List of all sujects IDs
Nsub = length(subjects);    % 144

cond_sets_new = {[1 3 5 7], [2 4 6 8]}; % {novid, vid}
set_name      = {'novid', 'vid'};

fband = [2 5];  % Hz

Nvox_full = 16102;  % full bilateral source-space size

% --- Get 'sensors' from first subject:
first_subject = cfg.subjects.include{1};
meg_pattern = fullfile(meg_dir, first_subject, '*_tsss_mc.fif');
fif_files = dir(meg_pattern);
if ~isempty(fif_files)
    meg_file = fullfile(meg_dir, first_subject, fif_files(1).name);
    raw = fiff_setup_read_raw(meg_file);
    sensors = get_sensors(raw);
end
%%
% 1-2. Per Subject: pool conditions, compute coherence (observed + surrogate)
%
% Pooling epochs across conditions BEFORE estimating the CSD/Coherence
% (rather than averaging four separately-estimated per-condition coherence
% maps) is the statistically correct order: coherence is a biased estimator
% whose bias depends on the number of epochs feeding the cross-spectral
% density, so pooling epochs irst gives one less biased estimate per
% condition set instead of averaging four differently-biased ones. This
% mirrors CM_combine_results uage for conditions 9 & 10, generalized here
% to conditions sets (1,3,5,7) and (2,4,6,8).

all_coh      = zeros(Nvox_full, 2, Nsub);   % observed
all_coh_perm = zeros(Nvox_full, 2, Nsub);   % surrogate
nave_used    = nan(Nsub, 2);
subjects_ok   = false(Nsub, 1);

for n_sub = 1:Nsub

    matfold      = fullfile(exp_dir, 'matfiles', subjects{n_sub});
    matfile_obs  = fullfile(matfold, [subjects{n_sub} '_main_coh.mat']);
    matfile_perm = fullfile(matfold, [subjects{n_sub} '_main_coh_perm_stat.mat']);

    if ~exist(matfile_obs, 'file') || ~exist(matfile_perm, 'file')
        warning('Missing coherence matfile for subject %s, skipping', subjects{n_sub});
        continue
    end

    S_obs  = load(matfile_obs, 'CMall');
    S_perm = load(matfile_perm, 'CMall');

    % --- forward model / inverse operator ---
    sub_fold_name = dir(fullfile(subjects_dir, ['meg', subjects{n_sub}(4:end) '*']));
    if isempty(sub_fold_name)
        warning('No freesurfer folder for %s, skipping', subjects{n_sub});
        continue
    end
    exp_name = 'ST';
    if ~strcmp(subjects{n_sub}(1:3), 'meg')
        exp_name = ['ST' subjects{n_sub}(1:3)];
    end
    fwdfile = fullfile(exp_dir, 'subjects', sub_fold_name(1).name, 'meg', ...
        [sub_fold_name(1).name '_from_MNI-' exp_name '-5-src-fwd.fif']);
    if ~exist(fwdfile, 'file')
        warning('No forward solution for %s, skipping', subjects{n_sub});
        continue
    end

    fwd = mne_read_forward_solution(fwdfile);
    L = zeros(length(sensors.picksgrads), fwd.sol.ncol/3*2);
    for n_source = 1:fwd.sol.ncol/3
        [U, ~, ~] = svd(fwd.sol.data(sensors.picksgrads, 3*n_source+(-2:0)), 'econ');
        L(:, 2*n_source+(-1:0)) = U(:, 1:2);
    end
    fwd.sol.data = L;
    fwd.sol.ncol = size(L, 2);
    fwd.sol.nrow = size(L, 1);
    fwd.source_nn(1:3:end, :) = [];
    fwd.nchan = length(sensors.picksgrads);

    COVnoise = sum(cat(3, S_obs.CMall.CSD), 3);
    COVnoise = COVnoise(sensors.picksgrads, sensors.picksgrads);
    CM_inv = MNE_inverse_MEEG(fwd, COVnoise);
    CM_inv.fwdfile = fwdfile;

    vertno_file = [fwdfile(1:end-4) '_vertno.mat'];
    if exist(vertno_file, 'file')
        load(vertno_file, 'vertno');
    else
        vertno = fwd.src.vertno; %#ok<NASGU>
        save(vertno_file, 'vertno');
    end

    Nvox = double(CM_inv.nsource);

    for n_set = 1:2
        conds_this_set = cond_sets_new{n_set};

        % --- pool epochs across conditions before computing coherence ---
        CM_obs  = S_obs.CMall(conds_this_set(1));
        CM_perm = S_perm.CMall(conds_this_set(1));
        for c = conds_this_set(2:end)
            CM_obs = CM_combine_results(CM_obs, S_obs.CMall(c));
            CM_perm = CM_combine_results(CM_perm, S_perm.CMall(c));
        end
        nave_used(n_sub, n_set) = CM_obs.nave;

        for which_data = 1:2
            if which_data == 1
                CM = CM_obs;
            else
                CM = CM_perm;
            end

            Find = find(CM.f >= fband(1) & CM.f <= fband(2));
            if isempty(Find)
                error('No frequency bins found in [%g;%g] Hz fos %s - check CM.freq_axis', fband(1), fband(2), subjects{n_sub});
            end

            Ngrid = 1000;
            teta = (1:Ngrid)'/Ngrid*pi;
            vec = [cos(teta) sin(teta)];

            coh = zeros(1, Nvox);
            for n_f = Find(:)'
                Cyy   = CM.Fyy(1, n_f, n_ref);
                T_CSD = CM_inv.Tinv * CM.CSD(sensors.picksgrads, sensors.picksgrads, n_f);
                Cs    = CM_inv.Tinv * CM.Fxy(sensors.picksgrads, n_f, n_ref);
                for n = 1:Nvox
                    Cxx = T_CSD([-1 0]+2*n, :) * CM_inv.Tinv([-1 0]+2*n, :)';
                    [Uc, Sd, Vc] = svd(Cxx);
                    if Sd(1,1)/Sd(2,2) > 100
                        Sd(2,2) = Sd(1,1)/100;  % regularize ill-conditionned direction
                        Cxx = Uc*Sd*Vc';
                    end
                    Cxy   = Cs([-1 0]+2*n);
                    num   = abs(vec*Cxy).^2;
                    denom = sum(vec*Cxx.*vec, 2)*Cyy;
                    coh(n) = coh(n) + max(real(num./denom));
                end
            end
            coh = coh / length(Find);

            if which_data == 1
                all_coh(vertno, n_set, n_sub) = coh;
            else
                all_coh_perm(vertno, n_set, n_sub) = coh;
            end
        end
    end
    subjects_ok(n_sub) = true;
    fprintf('Subject %s (%d/%d) done.\n', n_sub, Nsub);
end
save(fullfile(exp_dir, 'matfiles', 'all_coh_2-5Hz_mouth.mat'), 'all_coh', 'all_coh_perm', 'nave_used', 'subjects_ok', 'subjects', 'set_name', 'fband', 'n_ref');

                
%% 2.b Check epoch counts and whether nave correlates with age
%
% Coherence bias depends on the number of the epochs (nave) feeding the
% CSD.
% Pooling conditions changes per-subject nave relative to the original
% 10-conditions design, so this needs re-checking here rather than assuming
% previous code earlier nave-age check still appies.
load(fullfile(exp_dir, 'matfiles', 'all_coh_2-5Hz_mouth.mat'));
fprintf('\n --- Epoch counts (nave) summary per conditions set ---\n');
for n_set = 1:2
    fprintf("%s:    min: %d | median:%d | max:%d, N subjects:%d\n", set_name{n_set}, ...
    min(nave_used(subjects_ok, n_set)), median(nave_used(subjects_ok, n_set)), ...
    max(nave_used(subjects_ok, n_set)), sum(subjects_ok));
end
% Set an exclusion threshold - adjust based on the distribution above. As a
% starting point, exclude subjects whose pooled nave is markedly low (e.g.
% < half the median), since low nave inflates coherence bias and could act
% as a leverage point in the age correlation.
nave_threshold = 0.5 * median(nave_used(subjects_ok, :), 'all');
low_nave = any(nave_used < nave_threshold, 2) & subjects_ok;
if any(low_nave)
    fprintf('Flagging %d subject(s) with pooled nave below %.0f in at least one set.', sum(low_nave), nave_threshold);
    disp(subjects(low_nave));
else
    fprintf("No discrepency in pooled nave\n")
end
% Decide here whether to exclude flagged subjects before continuing.
% or carry them forward and treat as a sensitivity check.
include_sub = subjects_ok & ~low_nave;

% % Age vector aligned with 'subjects'
data_behav = '/mnt/usb-HardDrive_MB/expe_SpeechTrack/analysis/data_behav_dev.mat';
load (data_behav, 'behav_name', 'data_behav');
age = data_behav(:, 1);

for n_set = 1:2
    [r_nave, p_nave] = corr(age(include_sub), nave_used(include_sub, n_set), "type","Pearson");
    fprintf('%s: nave vs. age, Spearman r=%.2f, p=%.2f\n', set_name{n_set}, r_nave, p_nave);
end
% If either p is small, nave is a genuineconfound for this pooled design
% and should be reported for (e.g. partial correlation, or nave as a
% co-varaiate) rather than assumed neglible.

%% 4a. Grand average coherence maps -> localize occipital ROI (data driven)
%
% ROI is defined from the OBSERVED grand average only (never from the
% surrogate, and never from the age), so ROI selection is orthogonal to the
% age question - but the peak location is still chosen from the data, so
% absolute coherence values at that location will be somewhat inflated by
% chance (winner's curse). This is exactly why the subsequent age analysis
% uses the baseline-corrected value rather than the raw observed value: it
% does not remove the ROI-selection bias from the absolute coherence level,
% but it does prevent that bas (and any age-related SNR difference) from 
% leaking into the age correlation, since the same selection is applied
% identically to the observed and surrogates maps.

grand_avg = squeeze(mean(all_coh(: ,: ,include_sub), 3));   % Nvox x 2
% --- Sanity check: write grand_avg to NIfTI and inspevt in MRIcron ---
% Goal confirm the (still-to-be-found) peak sits inside a spatially
% extended, smooth-looking blob of elevated coherence spanning several
% neighboring vertices - not a single bright voxel surrounded by near-zero
% values, which would indicates noise rather than signal.
qc_fold = fullfile(exp_dir, 'subjects', 'source_maps_coh_2_5Hz_mouth', 'grand_avg_QC');
if ~exist(qc_fold, 'dir'); mkdir(qc_fold); end
mni_file = fullfile(exp_dir, 'subjects', 'MNI', 'bem', 'MNI-volume-5-src.fif');
for n_set = 1:2
    param = [];
    param.exp          = ['grand_avg_' set_name{n_set}];
    param.fold_fs      = 'ST';
    param.subjects_dir = fullfile(exp_dir, 'subjects');
    param.dxyz         = 1;
    param.w            = 1;
    param.srcfile      = mni_file;
    [map_mri, S]       = CM_overlay_MNI(grand_avg(:, n_set)', param);
    S.fname = fullfile(qc_fold, [param.exp '.nii']);
    S = spm_write_vol(S, map_mri);
end
%% 4.b Localize occipital maxim (ROI)
mni_file = fullfile(exp_dir, 'subjects', 'MNI', 'bem', 'MNI-volume-5-src.fif');
src_mni = mne_read_source_spaces(mni_file);
fprintf('src_mni.np = %d | Nvox_full expected: %d\n', src_mni.np, Nvox_full);

if src_mni.np == Nvox_full
    % rr already matches the used source indeing 1:1
    rr_mni = 1000 * src_mni.rr - [1 17 -19];
else
    % rr is a larger candidate grid - reuse the persubject vertno saved
    % earlier (identical across subjects since derived from the same
    % template forward computation) to subset it.
    any_sub = find(include_sub, 1);
    sub_fold_name = dir(fullfile(subjects_dir, ['meg' subjects{any_sub}(4:end) '*']));
    exp_name_chk = 'ST';
    if ~strcmp(subjects{any_sub}(1:3), 'meg')
        exp_name_chk = ['ST' subjects{any_sub}(1:3)];
    end

    fwdfile_chk = fullfile(exp_dir, 'subjects', sub_fold_name(1).name, 'meg', ...
        [sub_fold_name(1).name '_from_MNI-' exp_name_chk '-5-src-fwd.fif']);
    vertno_file_chk = [fwdfile_chk(1:end-4) '_vertno.mat'];
    load(vertno_file_chk, 'vertno');
    fprintf('length(vertno) = %d\n', length(vertno));
    rr_mni = 1000 * src_mni.rr(vertno, :);
end

% Occipital search mask: broad posterior constraint only, so the peak is
% still found by the data within that region rather than fixed by hand.
% Adjust the y-cutoff if it looks too permissive/restrictive once you see
% the grand-average map.
% prendre coord MNI de stat -> C'est mon nouveau mask

occ_mask = rr_mni(:,2) < -70 & rr_mni(:, 1) < 0;

roi_radius_mm = 15; % Cluster around the peak
ROI_vertices  = cell(1, 2);
peak_coord    = zeros(2, 3);
for n_set = 1:2
    map = grand_avg(:, n_set);
    map(~occ_mask) = -Inf;
    [peak_val, peak_idx] = max(map);
    peak_coord(n_set, :) = rr_mni(peak_idx, :);
    d = sqrt(sum((rr_mni - peak_coord(n_set, :)).^2, 2));
    ROI_vertices{n_set} = find(d <= roi_radius_mm);
    fprintf("Set %s: occipital peak coherence=%.3f at MNI [%.0f %.0f %.0f], ROI = %d vertices\n", ...
        set_name{n_set}, peak_val, peak_coord(n_set, 1), peak_coord(n_set, 2), peak_coord(n_set, 3), ...
        length(ROI_vertices{n_set}));
end
% Sanity check before proceeding: plot grand_avg on the cortical surface
% (e.g. via your usual CM_overlay_MNI + freeview/tksurfer pipeline) and
% confirm the peak looks like plausible, spatially coherent occipital
% cluster rather than an isolated noisy vertex.

%% 5. Extract per-subject occipital coherence (observed + surrogate)

roi_coh      = nan(Nsub, 2);
roi_coh_perm = nan(Nsub, 2);
for n_set = 1:2
    roi_coh(include_sub, n_set)      = squeeze(mean(all_coh(ROI_vertices{n_set}, n_set, include_sub), 1));
    roi_coh_perm(include_sub, n_set) = squeeze(mean(all_coh_perm(ROI_vertices{n_set}, n_set, include_sub), 1));
end

% Baseline-corrected value used for the used of age correlation
roi_coh_corrected = roi_coh ;%- roi_coh_perm; % à ne pas faire -> ROI coh 

save(fullfile(exp_dir, 'matfiles', 'roi_coh_2-5Hz_mouth_right.mat'), ...
    "roi_coh", "roi_coh_perm", "roi_coh_corrected", "ROI_vertices", "peak_coord", "subjects", "include_sub")

%% 6. Group-level statistical maps: observed vs.surrogate(per condition set)
%
% This establishes whether there is real, above-chance mouth-opening
% coherence anywhere in the brain (and specifically in the occipital ROI),
% which is the claim the ROI selection in step 4 depends on being true.
% Paired design because each subject contributes one observed and one
% surrogate map from the same epochs.

group_fold = fullfile(exp_dir, 'subjects', 'source_maps_coh_2-5Hz_mouth');
if ~exist(group_fold, 'dir'); mkdir(group_fold); end
 
stat = struct();
for n_set = 1:2
    sub_group_fold = fullfile(group_fold, set_name{n_set});
    perm_fold       = fullfile(sub_group_fold, 'perm_stat');
    if ~exist(sub_group_fold, 'dir'); mkdir(sub_group_fold); end
    if ~exist(perm_fold, 'dir'); mkdir(perm_fold); end
 
    for n_sub = find(include_sub)'
        param = [];
        param.exp         = ['w' subjects{n_sub}];
        param.fold_fs      = 'ST';
        param.subjects_dir = fullfile(exp_dir, 'subjects');
        param.dxyz         = 1;
        param.w            = 1;
        param.srcfile      = mni_file;
 
        [map_mri, S] = CM_overlay_MNI(all_coh(:, n_set, n_sub)', param);
        S.fname = fullfile(sub_group_fold, [param.exp '.nii']);
        S = spm_write_vol(S, map_mri);
 
        [map_mri_perm, Sp] = CM_overlay_MNI(all_coh_perm(:, n_set, n_sub)', param);
        Sp.fname = fullfile(perm_fold, [param.exp '.nii']);
        Sp = spm_write_vol(Sp, map_mri_perm);
    end
 
    cfg = [];
    cfg.niifold1   = sub_group_fold;
    cfg.niifiles1  = dir(fullfile(sub_group_fold, 'w*.nii'));
    cfg.niifold2   = perm_fold;
    cfg.niifiles2  = dir(fullfile(perm_fold, 'w*.nii'));
    cfg.filelabel1 = [set_name{n_set} '_observed'];
    cfg.filelabel2 = [set_name{n_set} '_surrogate'];
    cfg.smooth     = 8;
    cfg.fun        = @(x) x;
    cfg.alpha      = 0.05;
    cfg.design     = 'within';   % paired: each subject vs. own surrogate
 
    stat.(set_name{n_set}) = CM_stat_compare_source_maps_between(cfg);
end
save(fullfile(group_fold, 'stat_observed_vs_surrogate.mat'), 'stat');

%% 7. Correlate baseline-corrected occipital coherence with age
% 
% Spearman by default: robust to non-linearity/outliers, appropriate for a
% developmental sample without first assuming a linear relationship. Switch
% to Pearson only if the scatterplot below looks linear and homeoscedastic.
%
% Two tests here (novid, vid) - report which is primary, or apply a
% Holm/Bonferroni correction (alpha_corrected = 0.05/2) if both are
% presented as confirmatory rather than exploratory.

corr_results = struct();
for n_set = 1:2
    x = age(include_sub);
    y = roi_coh_corrected(include_sub, n_set);

    [r_s, p_s] = corr(x, y, 'type', 'Spearman');
    [r_p, p_p] = corr(x, y, 'type', 'Pearson');
    corr_results.(set_name{n_set}).spearman = [r_s, p_s];
    corr_results.(set_name{n_set}).pearson  = [r_p, p_p];

    figure;
    scatter(x, y, 40, 'filled'); lsline;
    xlabel('Age');
    ylabel(sprintf('Occipital coherence, right hemisphere (%s, 2-5 Hz mouth ref.)', set_name{n_set}));
    title(sprintf('%s: Spearman r=%.2f, p=%.3f | Pearson: r=%.2f, p=%.3f', ...
        set_name{n_set}, r_s, p_s, r_p, p_p));
end

disp(corr_results)

%% 8. Leave-one-out (non-circular) ROI redefinition and re-test
%
% WHY THIS STEP: the ROI in step 4 was located using the grand average
% across all included subjects - the same subjects whose values are then
% correlated with age. Each subject's own data contributed to placing the
% ROI, so a subject with an unusually high (or low) coherence value near a
% candidate peak can pull the peak toward their own dat point. This is
% "double-dipping" circular analysis (Kriegeskorte et al. 2009 "on voodoo
% correlation"): using the same data both to select a spatial feature and
% to test a hypothesis about it inflates the apparent effect and can
% produce a spurious "significant" correlation even if there is no true
% effect.
%
% Fix: leave-one-out (jacknife). For subject i, define the ROI using the
% grand average of everyone EXCEPT i, then extract i's value from that
% independently-defined location. Repeat for every subject. No subject's
% own data can then have influenced where "their" ROI sits, no any
% correlation that survives is not an artifact of the selection procedure.
% As a diagnostic, we also trackhow much the peak location moves across
% leave-one-out iterations: large jitter means the "peak" is unstable/noise
% dominated regardless of the correlation result; small jitter means a
% spatially consistent effect.

sub_idx = find(include_sub);
Nsub_incl = length(sub_idx);

roi_coh_loo      = nan(Nsub, 2);
roi_coh_perm_loo = nan(Nsub, 2);
peak_coord_loo   = nan(Nsub, 3, 2);

for n_set = 1:2
    for k = 1:Nsub_incl
        this_sub = sub_idx(k);
        others = sub_idx(sub_idx ~= this_sub);

        grand_avg_loo = mean(all_coh(:, n_set, others), 3);
        map = grand_avg_loo;
        map(~occ_mask) = -Inf;
        [~, peak_idx] = max(map);
        this_peak = rr_mni(peak_idx, :);
        peak_coord_loo(this_sub, :, n_set) = this_peak;

        d = sqrt(sum((rr_mni - this_peak).^2, 2));
        roi_loo = find(d <= roi_radius_mm);

        roi_coh_loo(this_sub, n_set)      = mean(all_coh(roi_loo, n_set, this_sub), 1);
        roi_coh_perm_loo(this_sub, n_set) = mean(all_coh_perm(roi_loo, n_set, this_sub), 1);
    end
end

roi_coh_corrected_loo = roi_coh_loo - roi_coh_perm_loo;

fprintf("\n--- ROI stability across leave-one-out iterations ---\n");
for n_set = 1:2
    d_from_full = sqrt(sum((peak_coord_loo(sub_idx, :, n_set) - peak_coord(n_set, :)).^2, 2));
    fprintf("%s: LOO peak displacement from full-sample peak: median=%.1fmm, max=%.1fmm\n", ...
        set_name{n_set}, median(d_from_full), max(d_from_full));
end
% Rule of thumb: displacement well within roi_radius_mm across most
% subjects suggests a stable location; displacements of several times
% roi_radius_mm for many subjects the "peak" is not robust spatial feature
% and any correlation should be treated cautiously regardless of its
% p-value.

fprintf("\n--- Leave-one-out (non-circular) age correlation ---\n");
corr_results_loo = struct();
for n_set = 1:2
    x = age(include_sub);
    y = roi_coh_corrected_loo(include_sub, n_set);

    [r_s, p_s] = corr(x, y, "type", "Spearman");
    [r_p, p_p] = corr(x, y, "type", "Pearson");
    corr_results_loo.(set_name{n_set}).spearman = [r_s, p_s];
    corr_results_loo.(set_name{n_set}).pearson = [r_p, p_p];

    figure;
    scatter(x, y, 40, 'filled'); lsline;
    xlabel('Age');
    ylabel(sprintf('LOO occipital coherence, baseline-corrected (%s, 2-5 Hz mouth ref.)', set_name{n_set}));
    title(sprintf('%s: Spearman r=%.2f, p=%.3f | Pearson: r=%.2f, p=%.3f', ...
        set_name{n_set}, r_s, p_s, r_p, p_p));

end

%% 9. Permutation (label-shuffle) test on the age correlation
%
% WHY THIS STEP: the parametric p-values above rely on assumptions -
% Pearsons assumes approximately bivariate-normal data, and Spearman's
% usual p_values is itself an asymptotic apporoximation that is less
% reliable as modest sample sizes, with ties, or with skewed distributions
% (coherence values are non-negative and often right-skewed). Your observed
% p-values sit right at the conventional significance edge, which is
% exactly where small deviations from these assumptions can flip a
% conclusion. Alabel-permutation test instead builds an empirical null
% distribution directly from your own data: repeatedly break the true
% subject-to-age pairing at random, recompute the correlation falls
% relative to that null. This requires no distributional assumption beyond
% exchangeability of subjects undr the null hypothesis of "no true
% association" - precisely the null you wan to rule out. We test the LOO
% (non-circular) values here, since those are ones you'd actually want to
% defend.

rng(1);
Nperm = 1000;
perm_test_results = struct();
corr_type = 'Spearman';
for n_set = 1:2
    x = age(include_sub);
    y = roi_coh_corrected_loo(include_sub, n_set);
    n = length(x);

    r_obs = corr(x, y, 'type', corr_type);

    r_null = nan(Nperm, 1);
    for p_i = 1:Nperm
        perm_idx = randperm(n);
        r_null(p_i) = corr(x(perm_idx), y, 'type', corr_type);
    end

    p_perm = mean(abs(r_null) >= abs(r_obs));   % Two sided empirical p

    perm_test_results.(set_name{n_set}).r_obs = r_obs;
    perm_test_results.(set_name{n_set}).p_perm = p_perm;


        fprintf('%s: observed %s r=%.3f, permutation p=%.4f (Nperm=%d)\n', ...
        set_name{n_set}, corr_type,r_obs, p_perm, Nperm);
 
    figure;
    histogram(r_null, 50); hold on;
    xline(r_obs, 'r', 'LineWidth', 2);
    xlabel(sprintf('Null %s r (age labels shuffled)', corr_type)); ylabel('Count');
    title(sprintf('%s: permutation null, observed r=%.3f, p=%.4f', set_name{n_set}, r_obs, p_perm));

end

% Multiple-comparisons correction across the two condition sets
% (Holm-Bonferroni on the two permutation p-values), consistent with the
% two-tests issue flagged since step 7.
p_vals = [perm_test_results.(set_name{1}).p_perm, perm_test_results.(set_name{2}).p_perm];
[p_sorted, order] = sort(p_vals);
holm_alpha = 0.05 ./ (2 - (0:1));
holm_reject = p_sorted < holm_alpha;
fprintf('\n--- Holm-Bonferroni correction across vid/novid permutation p-values ---\n');
for i = 1:2
        fprintf('%s: p=%.4f, threshold=%.4f, survives=%d\n', ...
        set_name{order(i)}, p_sorted(i), holm_alpha(i), holm_reject(i));
end

save(fullfile(exp_dir, 'matfiles', 'age_corr_robustness_checks.mat'), ...
    'roi_coh_corrected_loo', 'peak_coord_loo', 'corr_results_loo', 'perm_test_results', 'p_vals', 'holm_reject');





