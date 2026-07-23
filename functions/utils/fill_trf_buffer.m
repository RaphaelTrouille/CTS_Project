function trf_buf = fill_trf_buffer(trf_buf, cfg, bad, vid_en_in_SiN, t, dec, MISCorig, CM, raw, subj_files)
% FILL_TRF_BUFFER  Assign cross-val fold indices and accumulate trial data into trf_buf.
%
% Reproduces the randperm logic from config_speechtrack:
%   - For each of the n_folds-1 condition indices, find which temporal windows match
%   - If multiple windows match, draw one randomly (randperm)
%   - Last fold reuses the 2nd draw from fold n_folds-1's randperm
%
% Inputs:
%   trf_buf      - current buffer struct (accumulated across trials)
%   cfg          - pipeline config (cfg.trf.n_folds, cfg.conditions)
%   bad          - bad sample mask for this trial [1 x n_samples]
%   vid_en_in_SiN- condition target matrix used for window matching
%   t            - time vector of condition boundaries
%   dec          - sample offset (decimation/start index)
%   MISCorig     - struct with field Fs (sampling frequency)
%   CM           - struct with field ref (reference signals)
%   raw          - raw MEG object (used to extract sensor labels)
%   subj_files   - struct with field meg_file
%
% Output:
%   trf_buf      - updated buffer struct

this_cond_ind = zeros(size(bad));
n_t_sav       = [];
this_order    = [];

for n_fold = 1:cfg.trf.n_folds
    if n_fold == cfg.trf.n_folds
        % Reuse 2nd pick from previous fold's randperm
        n_t = n_t_sav(this_order(2));
    else
        n_t        = find(~any(bsxfun(@minus, cfg.conditions(n_fold).target, vid_en_in_SiN)));
        this_order = randperm(length(n_t));
        n_t_sav    = n_t;
        n_t        = n_t_sav(this_order(1));
    end

    t_start = max(dec + round( t(n_t)   * MISCorig.Fs + 1), 1);
    t_end   = min(dec + round( t(n_t+1) * MISCorig.Fs),    length(bad));
    this_cond_ind(t_start:t_end) = n_fold;
end

% Accumulate into buffer
trf_buf.meg_files{end+1} = subj_files.meg_file;
trf_buf.bad              = [trf_buf.bad,  bad];
trf_buf.cond_ind         = [trf_buf.cond_ind, this_cond_ind];
trf_buf.refs{end+1}      = CM.ref;
trf_buf.sensors          = get_sensors(raw);   % overwritten each trial (same sensors across trials)
trf_buf.Fs               = MISCorig.Fs;

end