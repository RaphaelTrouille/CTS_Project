function CMfwd = run_TRF_sensor(trf_buf, cfg)
% RUN_TRF_SENSOR  Forward TRF mapping in sensor space (MEG → audio envelope).
%
% DESCRIPTION:
%   Computes forward TRF models (evoked-like responses) in sensor space
%   using a 10-fold cross-validation scheme. For each frequency band
%   (n_type) and each fold (n_set), the model is trained on 9 folds and
%   tested on the held-out fold.
%
%   MEG data is read directly from .fif files (stored in trf_buf.meg_files)
%   rather than kept in RAM. Reference signals (audio envelopes) are
%   reconstructed from trf_buf.refs (CM.ref structs accumulated during
%   the trial loop in main_pipeline). Permutation is already applied to
%   CM.ref in main_pipeline when perm_stat is true — no shuffling here.
%
% INPUTS:
%   trf_buf   - Buffer struct accumulated during the trial loop. Fields:
%                 .meg_files  : cell array of .fif file paths (1 x Nvids)
%                 .bad        : concatenated artefact mask (1 x Nsamples)
%                 .cond_ind   : concatenated condition index (1 x Nsamples)
%                               0 = unassigned, 1-10 = cross-val fold index
%                 .refs       : cell array of CM.ref structs (1 x Nvids)
%                               each entry is a 1xNref struct with .chan and .info
%                               already shuffled if perm_stat was true in main_pipeline
%                 .sensors    : sensors struct from get_sensors (last trial)
%                 .Fs         : MEG sampling rate (Hz)
%   cfg       - Pipeline config struct. Relevant fields:
%                 .trf.freq_bands     : struct array with fields:
%                                         .label   : e.g. 'phrasal', 'syllabic'
%                                         .band    : [flo fhi] Hz
%                                         .width   : [flo fhi] transition Hz
%                                         .ds      : downsample factor
%                                         .t_pre   : pre-stimulus window (ms)
%                                         .t_post  : post-stimulus window (ms)
%                                         .t_buff  : buffer (ms)
%                                         .lambda  : regularisation vector
%                 .trf.refs_keep      : indices of CM.ref to use (e.g. [1 2 3 4])
%                 .trf.n_folds        : number of cross-val folds (default 10)
%                 .trf.att_ref_index  : index of attended envelope in CM.ref (default 2)
%                                       used to compute speak_time silence mask
%
% OUTPUTS:
%   CMfwd - Forward TRF results, shape (n_folds x n_types).
%           Each entry is a CM struct after CM_TRF_MEEG_train with fields
%           including .w (weights), .t (time axis), .lambda, .n_l, etc.
%
% NOTES:
%   - speak_time (silence mask for syllabic band) is computed here from
%     the attended envelope but is reserved for the backward mapping —
%     not used in the forward pass.
%   - Forward mapping runs on all gradiometers (picks_sens{3}).
%   - Left/right hemisphere splits are reserved for backward mapping.
%   - n_pass == 2 (TRF time course extraction) not yet implemented.
%
% DEPENDENCIES:
%   - fiff_setup_read_raw        (MNE/FieldTrip)
%   - fiff_read_raw_segment      (MNE/FieldTrip)
%   - CM_TRF_MEEG_train          (cartographie_motrice toolbox)
%   - log_msg                    (custom)
%
% USAGE:
%   CMfwd = run_TRF_sensor(trf_buf, cfg);
% =========================================================================

Fs       = trf_buf.Fs;
sensors  = trf_buf.sensors;
Nvids    = length(trf_buf.meg_files);
n_folds  = cfg.trf.n_folds;

%% 1. LOAD MEG DATA & RECONSTRUCT ENVELOPES
% Read gradiometer data from each .fif file and stack
log_msg(cfg, '  [TRF] Loading MEG data from %d trials...\n', Nvids);

data = [];

for n_vid = 1:Nvids
    raw       = fiff_setup_read_raw(trf_buf.meg_files{n_vid});
    this_data = fiff_read_raw_segment(raw, 0, inf, ...
                    sensors.picksMEG(sensors.picksgrads));
    data      = cat(2, data, this_data);
end

%% 2. TEMPORAL CONCATENATION OF REFERENCE SIGNALS
% Concatenates trial-based references (from trf_buf) into a single continuous 
% matrix (Yenv = Refs x t).
Nrefs    = length(trf_buf.refs{1});
Nsamples = size(data, 2);
Yenv     = zeros(Nrefs,Nsamples);
offset   = 0;
for n_vid = 1:Nvids
    refs   = trf_buf.refs{n_vid};
    n_samp = size(refs(1).chan, 2);
    for n_ref = 1:Nrefs
        Yenv(n_ref, offset +(1:n_samp)) = refs(n_ref).chan;
    end
    offset = offset + n_samp;
end

bad      = trf_buf.bad;
cond_ind = trf_buf.cond_ind;

%% 3. SPEAK_TIME MASK
% Marks low energy (silence) periods as 0 based on attended speech envelope.
% Reserved for backward mapping - not used in forward pass.
att_idx            = cfg.trf.att_ref_index;
realenv            = Yenv(att_idx, :);
realenv(realenv == 0) = Nan;
thlow              = prctile(realenv, 5) * 5;
realenv(isnan(realenv)) = 0;
deb                = find(diff(realenv < thlow) == 1) + 1;
fin                = find(diff(realenv < thlow) == -1);
if fin(1) < deb(1), deb = [1 deb]; end
if deb(end) > fin(end), fin = [fin length(realenv)]; end
fin                = max(fin - 0.1*Fs, 1);
deb                = min(deb + 0.1*Fs, length(realenv));
to_rm              = find(fin - deb < 0.2*Fs);
deb(to_rm)         = [];
fin(to_rm)         = [];
speak_time         = ones(1, length(realenv));
for n_deb = 1:length(deb)
    speak_time(deb(n_deb):fin(n_deb)) = 0;
end

%% 4. FORWARD TRF - CROSS-VALIDATION LOOP

refs_keep = cfg.trf.refs_keep;
CMfwd     = [];

for n_type = 1:length(cfg.trf.freq_bands)
    band = cfg.trf.freq_bands(n_type);
    log_msg(cfg, '  [TRF] Band %d/%d: %s\n', n_type, length(cfg.trf.freq_bands), band.label);

    % Build CM template for this frequency band
    CMtemp              = [];
    CMtemp.map          = 1;
    CMtemp.Fs           = Fs;
    CMtemp.Nfold        = n_folds;
    CMtemp.filt.par     = {'high', 'low'};
    CMtemp.filt.f_vect  = band.band / Fs * 2;
    CMtemp.filt.Ws_vect = band.width / Fs * 2;
    CMtemp.ds           = band.ds;
    CMtemp.t_pre        = band.t_pre;
    CMtemp.t_post       = band.t_post;
    CMtemp.t_buff       = band.t_buff;
    CMtemp.lambda       = band.lambda;

    % Assign normalised reference signals to CM template
    for n_ref = 1:Nrefs
        sig = Yenv(n_ref, :);
        sig = sig - mean(sig);
        sig = sig / max(std(sig));
        CMtemp.ref(n_ref).chan = sig;
        CMtemp.ref(n_ref).info = trf_but.refs{1}(n_ref).info;
        CMtemp.ref(n_ref).filt = [];
        CMtemp.ref(n_ref).norm = 0;
    end

    for n_set = 1:n_folds
        
        % Training set: exclude current fold + unassigned (cond_in == 0)
        this_bad_train = bad;
        for n_cond = [0 n_set]
            this_bad_train(cond_ind == n_cond) = 1;
        end

        % Expand clean segments with 4s buffer on each side
        tmp = diff([1 this_bad_train 1]);
        tdebs = find(tmp == -1);
        tfins = find(tmp == 1) - 1;
        to_keep = zeros(size(bad));
        for n = 1:length(tdebs)
            to_keep(max(tdebs(n) - 4*Fs, 1) : min(tfins(n) + 4*Fs, length(bad))) = 1;
        end
        to_keep_train = find(to_keep);

        % Train
        this_CM = CMtemp;
        for n_ref = 1:length(CMtemp.ref)
            this_CM.ref(n_ref).chan = CMtemp.ref(n_ref).chan(:, to_keep_train);
        end
        this_CM.ref = this_CM.ref(n_ref).chan(:, to_keep_train);
        this_CM     = CM_TRF_MEEG_train(this_CM, ...
                                        this_bad_train(to_keep_train), ...
                                        data(1:size(data,1), to_keep_train));
        this_CM =rmfield(this_CM, {'ref' 'tdeb' 'tfin'});

        % Test TO DO

        
        if n_set == 1 && n_type == 1
            CMfwd = this_CM;
        else
            CMfwd(n_set, n_type) = this_CM;
        end

        log_msg(cfg, '  fold %d/%d done\n', n_set, n_folds);

    end  % n_set
end  % % n_type

end