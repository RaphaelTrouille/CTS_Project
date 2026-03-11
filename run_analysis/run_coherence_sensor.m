function CMall = run_coherence_sensor(CMall, CM, bad, n_vid, n_cond, cfg)
% RUN_COHERENCE_SENSOR  Run sensor-space cortical coherence for one
%                        condition and accumulate results into CMall.
%
% DESCRIPTION:
%   Wraps CM_coh_MEG_ref_one_pass for one trial/condition, handles optional
%   epoch-count equalisation across conditions, cleans up intermediate fields,
%   zeros out results when no valid epochs were found, and accumulates the
%   result into CMall using CM_combine_results.
%
% INPUTS:
%   CMall  - Accumulated results array (may be empty on first call).
%            Updated and returned.
%   CM     - CM structure for the current trial/condition, containing:
%              .infile, .ref, .tdeb, .tfin, .quantum, .win, .Find, etc.
%   bad    - Binary artefact mask (1 x N); 1 = artefact, 0 = clean
%   n_vid  - Current video/trial index (used for epoch equalisation logic)
%   n_cond - Current condition index
%   cfg    - Pipeline configuration struct (see config.m). Relevant fields:
%              .cm.equalize_nave : (bool) equalise epoch count across conditions
%              .cm.n_rm          : (1 x Ncond) epochs to remove per condition
%              .cm.equalize_vids : (1 x 2) [vid_index, max_cond] defining
%                                  which trials/conditions to equalise
%                                  e.g. [2 4; 3 4] means:
%                                    vid2: equalise if n_cond <= 4
%                                    vid3: equalise if n_cond > 4
%
% OUTPUTS:
%   CMall  - Updated results array (1 x Ncond), where each entry is the
%            accumulated coherence across trials for that condition
%
% NOTES:
%   - EEG coherence (CM_coh_EEG_ref_one_pass) is not implemented here.
%     Add a cfg.space flag and a separate branch if needed.
%   - Epoch equalisation runs CM_coh_MEG_ref_one_pass a second time with
%     CM.maxave set to limit the number of epochs, then removes the field.
%   - Fields 'fxy', 'spectre', and 'ref' are removed before accumulation
%     to keep memory usage low. CSD is zeroed if nave == 0.
%
% DEPENDENCIES:
%   - CM_coh_MEG_ref_one_pass  (cartographie_motrice toolbox)
%   - CM_combine_results       (cartographie_motrice toolbox)
%
% USAGE:
%   CMall = run_coherence_sensor(CMall, CM, bad, n_vid, n_cond, cfg);
%
% -------------------------------------------------------------------------

    %% 1. COMPUTE COHERENCE
    CMresult = CM_coh_MEG_ref_one_pass(CM, bad);
    
    %% 2. OPTIONAL EPOCH COUNT EQUALISATION
    % Re-run with a capped maxave to equalise the number of epochs cross
    % conditions. Useful when trial duration differs between conditions.
    % Controlled bu cfg.cm.equaliza_nave and cfg.cm.equalize_vids in config.m.
    if cfg.cm.equalize_nave && ~isempty(cfg.cm.n_rm)
        if should_equalize(n_vid, n_cond, cfg.equalize_vids)
            CM.maxave = CMresult.nave - cfg.cm.n_rm(n_cond);
            CMresult  = CM_coh_MEG_ref_one_pass(CM, bad);
            CM        = rmfield(CM, 'maxave'); % Clean up so it doesn't persist
        end
    end
    
    %% 3. CLEAN UP INTERMEDIATE FIELDS
    % Remove large fields not needed for group-level analysis
    fields_to_remove = {'fxy', 'spectre', 'ref'};
    for f = 1:length(fields_to_remove)
        if isfield(CMresult, fields_to_remove{f})
            CMresult = rmfield(CMresult, fields_to_remove{f});
        end
    end
    
    %% 4. ZERO OUT IF NO VALID EPOCHS
    if CMresult.nave == 0
        CMresult.Fxx(:)     = 0;
        CMresult.Fyy(:)     = 0;
        CMresult.Fxy(:)     = 0;
        CMresult.Fgrads(:)  = 0;
        if isfield(CMresult, 'CSD')
            CMresult.CSD(:) = 0;
        end
    end
    
    %% 5. ACCUMULATE INTO CMall
    if isempty(CMall) || size(CMall, 2) < n_cond
        % First trial for this condition - initialise
        CMall(n_cond) = CMresult;
    else
        % Subsequent trials - merge with existing results
        CMall(n_cond) = CM_combine_results(CMall(n_cond), CMresult);
    end
end


% -------------------------------------------------------------------------
% LOCAL HELPER
% -------------------------------------------------------------------------
function tf = should_equalize(n_vid, n_cond, equalize_vids)
% Check wether epoch equalisation should be applied for this vid/ond pair.
% equalize_vids is a M x 2 matrix: each row is [vid_index, threshold_cond]
%   - if threshold > 0 : equalie when n_cond <= threshold
%   - if threshold < 0 : equalise when n_cond > 0 abs(threshold)
%
% Example in config.m
%   cfg.cm.equalize_vids = [2,  4;  % vid2: equalise if n_cond <= 4
%                           3, -4]  % vid3: equalise if n_cond > 4

    tf = false;
    if isempty(equalize_vids), return; end

    for k = 1:size(equalize_vids, 1)
        vid_k = equalize_vids(k, 1);
        thr_k = equalize_vids(k, 2);
        if n_vid == vid_k
            if thr_k > 0 && n_cond <= thr_k
                tf = true; return
            elseif thr_k < 0 && n_cond > abs(thr_k)
                tf = true; return
        
            end
        end
    end
end