% =========================================================================
% MAIN_PIPELINE.M _ SpeechTrack CTS analysis pipeline
% =========================================================================
%
% USAGE: 
%   Edit config.m to set parameters, then run this script.
%
% PIPELINE STRUCTURE:
%   0. Environment setup & config loading
%   1. Subject loop
%      ├── Exclusion check
%      ├── Output path & skip logic
%      └── Trial loop (per video)
%           ├── A. Data loading & MEG/audio alignment
%           ├── B. Preprocessing (bandpass, reref, artifact rejection, zscore)
%           ├── C. Reference signal construction
%           └── D. Condition loop → analysis dispatch
%                ├── Coherence
%                ├── TRF
%                ├── ERP
%                └── Beamforming
%
% TO ADD A NEW ANALYSIS:
%   1. Add cfg.analysis.my_analysis = true/false in config.m
%   2. Add parameters in config.m
%   3. Implement run_my_analysis(CM, bad, cfg) at the bottom of this file
%   4. Add an elseif block in Section D below
%
% =========================================================================

%% 0. ENVIRONMENT & CONFIGURATION
% -------------------------------------------------------------------------
scriptPath = fileparts(mfilename("fullpath"));
if ~exist('setup_environment.m', 'file')
    addpath(fullfile(scriptPath, '..'));
end

[meg_dir, subjects, deriv_dir, snd_dir, vid_dir] = setup_environment();
config;     % Loads cfg struct - edit config.m to change all parameters

% Pre-compute MEL center frequencies (shared across all subjects/trials)
cf = get_mel_freqs(cfg.freq.mel_low, cfg.freq.mel_high, cfg.freq.mel_bands);

% Initialise results collectors
gof         = nan(length(subjects), cfg.expected_trials);
failed_subs = {};

%% 1. SUJECT LOOP
% =========================================================================
perm_passes = iif(cfg.perm_stat, [false true], false);

for perm_stat = perm_passes
    perm_label = iif(perm_stat, '[PERM]', '[MAIN]');

    for n_sub = 1:length(subjects)
        tic
        sub_name = subjects(n_sub).name;

        % -----------------------------------------------------------------
        % 1.1 EXCLUSION CHECK
        % -----------------------------------------------------------------
        [skip, skip_msg] = check_exclusion_criteria(sub_name,...
            cfg.subjects.exclude_ica, cfg.subjects.exclude_diagnosis);
        if skip
            log_msg(cfg, '%s [SKIP] %s - \n', perm_label, sub_name, skip_msg);
            continue
        end

        % -----------------------------------------------------------------
        % 1.2 OUTPUT PATH & SKIP LOGIC
        % -----------------------------------------------------------------
        
        mat_path = fullfile(deriv, 'results', sub_name);
        if ~exist(mat_path, 'dir'), mkdir(mat_path); end

        suffix    = iif(perm_stat, '_perm', '_main');
        save_file = fullfile(mat_path, [sub_name '_CTS_results' suffix '.mat']);

        if ~cfg.overwrite_results && exist(save_file, 'file')
            log_msg(cfg, '%s [SKIP] %s - result filealready exists.\n', perm_label, sub_name)
            continue
        end

        log_msg(cfg, '\n%s [%d/%d] Processing: %s\n', perm_label, n_sub, length(subjects), sub_name);

        % -----------------------------------------------------------------
        % 1.3 LOCATE SUBJECT NEURO DATA FOLDER (hendles dated subfolders)
        % -----------------------------------------------------------------

    end  % subjects loop
end  % perm loop