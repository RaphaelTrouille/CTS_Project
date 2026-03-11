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
%   2. └── Trial loop (per video)
%           ├── A. Data loading & MEG/audio alignment
%           ├── B. Reference signal construction
%           ├── C. MEG metadata & artefact mask
%           └── D. Condition loop → analysis dispatch
%                ├── Coherence
%                ├── TRF
%                └── ERP
%   3. Summary & Save the results
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
[meg_dir, subjects, deriv_dir, snd_dir, vid_dir] = setup_environment();
cfg = cts_config();     % Loads cfg struct - edit config.m to change all parameters

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
        mat_path  = fullfile(deriv_dir, 'results', sub_name);
        if ~exist(mat_path, 'dir'), mkdir(mat_path); end

        suffix    = iif(perm_stat, '_perm', '_main');
        save_file = fullfile(mat_path, [sub_name '_CTS_results' suffix '.mat']);

        if ~cfg.overwrite_results && exist(save_file, 'file')
            log_msg(cfg, '%s [SKIP] %s - result file already exists.\n', perm_label, sub_name)
            continue
        end

        log_msg(cfg, '\n%s [%d/%d] Processing: %s\n', perm_label, n_sub, length(subjects), sub_name);

        % -----------------------------------------------------------------
        % 1.3 LOCATE SUBJECT NEURO DATA FOLDER (hendles dated subfolders)
        % -----------------------------------------------------------------
        subfold_root = fullfile(meg_dir, sub_name);
        date         = dir(fullfile(subfold_root, '2*'));
        if isempty(date)
            subfold  = subfold_root;
        else
            subfold  = fullfile(subfold_root, date(1).name);
            log_msg(cfg, '  [INFO] Dated subfolder found, using root: %s. \n', subfold);
        end
        
        % -----------------------------------------------------------------
        % 1.4 GET STIMULUS SET ORDER & VIDEO LIST
        % -----------------------------------------------------------------
        try
            [set, order, vids] = get_set_order_vids(subfold, cfg.expected_trials);
        catch err
            log_msg(cfg, '  [WARNING] Could not retrieve set/order/vids: %s\n', err.message);
            failed_subs{end+1} = sub_name;
        end

        CMall = [];

        %% 2. TRIAL LOOP
        % =================================================================
        for n_vid = 1:length(vids)
            log_msg(cfg, '  > Trial %d/%d (vid %d)\n', n_vid, length(vids), vids(n_vid));

            subj_files = get_subject_files(subfold, snd_dir, sub_name, set, order, vids(n_vid));

            % Check that MEG file exists before going further
            if ~exist(subj_files.meg_file, 'file')
                log_msg(cfg, '  [WARNING] MEG file not found: %s - skipping trial. \n', subj_files.meg_file);
                continue
            end

            % Retrieve trial timings and condition flags
            [t, vid_en_in_SiN, t_dis] = get_vid_timings(order, vids(n_vid));

            %% A. DATA LOADING & ALIGNMENT
            % -------------------------------------------------------------
            Yglb     = load_WAV_audio(subj_files.snd_global);
            MISCorig = load_MISC(subj_files.meg_file, perm_stat);

            % Sanity check: audio duration must match theoretical timings
            if abs(t(end) - length(Yglb.signal) / Yglb.Fs) > 0.1
                log_msg(cfg, '  [WARNING] Audio/timing mismatch for vid%d - skipping trial. \n', vids(n_vid));
                continue
            end

            % Align MEG and audio (load cached result availbale)
            [dec, tds] = realign_sound_file(subj_files.sync_file, MISCorig, Yglb);
            
            % Verify and log alignment quality
            [t1, t2, L, gof_val] = alignement_verification(tds, dec, MISCorig, Yglb, perm_stat, false);
            gof(n_sub, n_vid) = gof_val;

            if ~perm_stat && ~isempty(gof_val) && gof_val < 0.5
                log_msg(cfg, '  [WARNING] Poor alignment (gof=%.2f) for vid%d.\n', gof_val, vids(n_vid));
            end

            %% B. REFERENCE SIGNAL CONTRUCTION
            % -------------------------------------------------------------
            CM = init_CM(MISCorig.Fs, cfg.cm.quantum, cfg.cm.window_type);

            for i = 1:size(cfg.ref_sources, 1)
                ref_label    = cfg.ref_sources{i, 1};
                ref_type     = cfg.ref_sources{i, 2};
                ref_file_key = cfg.ref_sources{i, 3};
            
                % Build envelope / motor signal
                switch ref_type
    
                    case 'audio'
                        wav = load_WAV_audio(subj_files.(ref_file_key));
                        Yp  = envelope_extraction(wav, cf);
    
                    case 'lips'
                        wav = load_WAV_audio(subj_files.snd_att);
                        Yp = lips_apperture(vid_dir, vids(n_vid), wav);

                    case 'envelope'
                        % Pre-computed low-frequency signal - add loading
                        % logic here
                        log_msg(cfg, '  [WARNING] Envelope source "%s" needs custom loading - skipping.\n', ref_label);
                        continue
                    
                    otherwise
                        log_msg(cfg, '  [WARNING] Unknown ref type "%s" - skipping.\n')
                        continue
                end
            

                % Align and downsample to MEG rate
                MISC = downsample_signal(Yp, tds, t1, t2, size(MISCorig.signal));
    
                % Permute if in surrogate mode
                if perm_stat
                    MISC = permute_signal_segments(MISC, MISCorig.fs);
                end
    
                CM = add_CM_ref(CM, ref_label, MISC);

            end % references loop
            
            clear Yglb MISC Yp wav % Free memory - no longer needed
            
            %% C. MEG METADATA & ARTEFACT MASK
            % -------------------------------------------------------------
            % % prepare_CM_data sets CM.infile, CM.Fs, CM.first_samp,
            % CM.last_samp — everything CM_coh_MEG_ref_one_pass needs
            % to read the raw data itself. It also builds the bad mask.

            [CM, bad] = prepare_CM_data(CM, dec, L, subj_files.meg_file);

            % Flag distractor/singing segments as artefacts
            bad = flag_distractor_segments(bad, t_dis, tds, dec, MISCorig.Fs);

            %% D. CONDITION LOOP -> ANALYSIS DISPATCH
            % -------------------------------------------------------------

            for n_cond = 1:length(cfg.conditions)
                cond   = cfg.conditions(n_cond);
                n_t    = find(~any(bsxfun(@minus, cond.target, vid_en_in_SiN)));
                
                if isempty(n_t)
                    log_msg(cfg, '  [INFO] Condition %d: %s not found in the trial: vid%d\n',...
                            n_cond, cond.label, vids(n_vid));
                    continue % This condition is not present in the current trial
                end

                log_msg(cfg, '  Condition %d: %s (%d windows)\n', n_cond, cond.label, length(n_t));
            
                % Set analysis time window (absolute sample indices)
                CM.tdeb = CM.first_samp + max(dec + t(n_t) * MISCorig.Fs, 0);
                CM.tfin = CM.first_samp + min(dec + t(n_t+1) * MISCorig.Fs, ...
                                               CM.last_samp - CM.first_samp);

                % ---------------------------------------------------------                                           CM.last_samp - CM.first_samp);
                % ANALYSIS DISPATCH
                % To add a new analysis: implement run_X(CM, bad, cfg) below
                % and add an elseif block here.
                % ---------------------------------------------------------
                
                CMresult = [];
    
                if cfg.analysis.coherence && cfg.space.sensor
                    CMall = run_coherence_sensor(CMall, CM, bad, n_vid, n_cond, cfg);

                elseif cfg.analysis.coherence && cfg.space.source
                    % TODO
                    continue

                elseif cfg.analysis.TRF && cfg.space.source
                    % TODO
                
                else
                    log_msg(cfg, '  [WARNING] No active analysis matched - check config.m\n');
                    continue
                end
    
                if isempty(CMall), continue; end
    
   
            end % conds loop
        end  % trials loop

        %% 3. SAVE RESULTS
        % -----------------------------------------------------------------
        % Extract surrogate statistics from CMall before saving.
            % CM_coh_MEG_ref_one_pass stores surrogate results in CMall.surrog
            % when CM.surrog is defined. Here we pull them out into a separate
            % array (surrog) indexed by [condition x trial], and clear them
            % from CMall to keep the saved file compact.
        if ~isempty(CMall)

            % Extract surrogate statistics from CMall before saving.
            % CM_coh_MEG_ref_one_pass stores surrogate results in CMall.surrog
            % when CM.surrog is defined. Here we pull them out into a separate
            % array (surrog) indexed by [condition x trial], and clear them
            % from CMall to keep the saved file compact.
            surrog = [];
            if isfield(CMall(1), 'surrog')
                for s1 = 1:size(CMall, 1)
                    for s2 = 1:size(CMall, 2)
                        surrog(s1, s2) = CMall(s1, s2).surrog;
                        CMall(s1, s2).surrog = [];
                    end
                end
            end

            gof_align = gof(n_sub, :);
            save(save_file, 'CMall', 'surrog', 'gof_align', 'cfg');
            log_msg(cfg, '  [SAVED] %s\n', save_file);
        else
            log_msg(cfg, '  [WARNING] CMall is empty — nothing saved for %s.\n', sub_name);
        end

        elapsed = toc;
        log_msg(cfg, '  Done in %.1f sec.\n', elapsed);

    end  % subjects loop
end  % perm loop

% Summary 
if ~isempty(failed_subs)
    fprintf('\n[SUMMARY] %d subject(s) failed:\n', length(failed_subs));
    fprintf('  - %s\n', failed_subs{:});
end
fprintf('\n[DONE] Pipeline complete.\n')