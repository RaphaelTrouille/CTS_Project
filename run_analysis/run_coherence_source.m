function run_coherence_source(CMall, n_cond, sub_name, perm_stat, cfg)
% RUN_COHERENCE_SOURCE  Project sensor-space coherence results to source
%                       space using DICS beamforming for one condition.
%
% DESCRIPTION:
%   Takes an already-computed CMall (output of run_coherence_sensor) and
%   projects the CSD matrix to source space using a subject-specific
%   forward solution. For each configured frequency band, builds a CM_CSD
%   structure, runs CM_DICS_MEEG to produce FreeSurfer overlay files (.w),
%   and moves them to the group results folder.
%
%   This function produces files on disk and has no return value.
%   It must be called after run_coherence_sensor has populated CMall.CSD.
%
% INPUTS:
%   CMall     - Accumulated sensor-space coherence results (1 x Ncond).
%               Must contain .CSD, .f, .Fs, .quantum, and .infile fields.
%               (requires CM.CSD = [] at init_CM time to trigger CSD computation)
%   n_cond    - Index of the condition to project (into CMall)
%   sub_name  - Subject identifier string (e.g. 'Meg3877')
%   perm_stat - Boolean; if true, appends '_perm_stat' to output filenames
%   cfg       - Pipeline configuration struct. Relevant fields:
%                 .beam.fwd_file     : forward solution path template.
%                                      Use {sub_name} as placeholder, e.g.:
%                                      '/subjects/{sub_name}/meg/{sub_name}_fwd.fif'
%                 .beam.subjects_dir : FreeSurfer subjects directory
%                 .beam.fold_fs      : FreeSurfer folder name for overlays
%                 .beam.dxyz         : source grid spacing (mm)
%                 .beam.freq_bands   : struct array with fields:
%                                        .label : band name (e.g. 'phrasal')
%                                        .Find  : frequency indices into CM.f
%                 .beam.ref_index    : reference signal index for CSD (default 2)
%                 .beam.group_fold   : root output folder for group overlays
%                 .beam.n_sensors    : number of MEG sensors (default 306)
%                 .conditions(n_cond).label : condition label for filenames
%
% OUTPUTS:
%   (none) — writes FreeSurfer overlay files (.w) to cfg.beam.group_fold
%
% NOTES:
%   - SVD reduction keeps the 2 dominant orientations per source dipole,
%     applied to gradiometer channels only.
%   - Noise covariance is estimated by averaging CSD across all conditions
%     and all frequencies >= 1 Hz.
%   - Output overlays are moved (not copied) from the FreeSurfer subject
%     folder to the group folder via unix mv.
%
% DEPENDENCIES:
%   - mne_read_forward_solution  (MNE/FieldTrip)
%   - MNE_inverse_MEEG           (cartographie_motrice toolbox)
%   - CM_DICS_MEEG               (cartographie_motrice toolbox)
%   - fiff_setup_read_raw        (MNE/FieldTrip)
%   - get_sensors                (custom)
%   - iif                        (custom)
%
% USAGE:
%   run_coherence_source(CMall, n_cond, sub_name, perm_stat, cfg);
%
% =========================================================================
    Nsens = cfg.beam.n_sensors;  % Number of MEG sensors (default 306)
     
    %% 1. SENSOR-SPACE COHERENCE (with CSD)
    % CM.CSD must be initialized to [] in init_CM for CSD to be computed
    CM_cond = CMall(n_cond);
     
    %% 2. NOISE COVARIANCE FROM CSD ACROSS ALL CONDITIONS
    % Average CSD across all conditions and frequencies >= 1 Hz as noise estimate
    COVnoise = real(cat(4, CMall.CSD));
    f_idx    = find(CM_cond.f >= 1);
    COVnoise = mean(mean(COVnoise(1:Nsens, 1:Nsens, f_idx, :), 3), 4);
     
    %% 3. FORWARD SOLUTION — SVD REDUCTION TO 2 ORIENTATIONS PER SOURCE
    % Build subject-specific forward solution path
    fwdfile = strrep(cfg.beam.fwd_file, '{sub_name}', sub_name);
    fwd     = mne_read_forward_solution(fwdfile);
     
    % Read sensor info to get gradiometer picks
    raw     = fiff_setup_read_raw(CM_cond.infile);
    sensors = get_sensors(raw);
     
    % Reduce each source dipole to its 2 dominant orientations via SVD
    % (applied to gradiometer channels only)
    L = zeros(length(sensors.picksgrads), fwd.sol.ncol / 3 * 2);
    for n_source = 1:double(fwd.sol.ncol) / 3
        [U, ~, ~] = svd(double(fwd.sol.data(sensors.picksgrads, 3*n_source+(-2:0))), 'econ');
        L(:, 2*n_source+(-1:0)) = U(:, 1:2);
    end
     
    % Update forward solution with reduced leadfield
    fwd.sol.data   = L;
    fwd.sol.ncol   = size(L, 2);
    fwd.sol.nrow   = size(L, 1);
    fwd.source_nn(1:3:end, :) = [];
    fwd.nchan      = length(sensors.picksgrads);
     
    %% 4. MNE INVERSE OPERATOR
    CM_inv          = MNE_inverse_MEEG(fwd, COVnoise(sensors.picksgrads, sensors.picksgrads));
    CM_inv.fwdfile  = fwdfile;
     
    %% 5. SOURCE PROJECTION PER FREQUENCY BAND
    for n_band = 1:length(cfg.beam.freq_bands)
        band  = cfg.beam.freq_bands(n_band);
        Find  = band.Find;
        ref   = cfg.beam.ref_index;  % Reference signal index (default 2 = attended speech)
     
        % Build CSD structure for this frequency band
        CM_CSD       = [];
        CM_CSD.Fs    = CM_cond.Fs;
        CM_CSD.F     = (Find - 1) * CM_CSD.Fs / CM_cond.quantum;
        CM_CSD.Find  = Find;
        CM_CSD.bads  = [];
        CM_CSD.CSD   = CM_cond.CSD([sensors.picksgrads, Nsens+ref], ...
                                   [sensors.picksgrads, Nsens+ref], Find);
     
        % Build output parameters
        param              = [];
        exp_label          = [sub_name '_' cfg.conditions(n_cond).label '_' band.label];
        param.exp          = iif(perm_stat, [exp_label '_perm_stat'], exp_label);
        param.fold_fs      = cfg.beam.fold_fs;
        param.subjects_dir = cfg.beam.subjects_dir;
        param.subject      = sub_name;
        param.dxyz         = cfg.beam.dxyz;
        param.w            = 1;  % Do not normalise back to MRI
     
        % Run DICS source projection
        CM_DICS_MEEG(param, CM_inv, CM_CSD);
     
        %% 6. MOVE OVERLAYS TO GROUP FOLDER
        sub_fold   = fullfile(cfg.beam.subjects_dir, sub_name, param.fold_fs);
        group_fold = fullfile(cfg.beam.group_fold, cfg.conditions(n_cond).label, band.label);
        if perm_stat
            group_fold = fullfile(group_fold, 'perm_stat');
        end
        if ~exist(group_fold, 'dir')
            mkdir(group_fold);
        end
        unix(['mv ' fullfile(sub_fold, ['w' param.exp '*']) ' ' group_fold]);
    end
end