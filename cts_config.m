function cfg = cts_config()
% CTS_CONFIG  Return the central configuration structure for the
%             SpeechTrack CTS pipeline.
%
% DESCRIPTION:
%   This is the ONLY file you need to edit between runs.
%   All analysis parameters, subject lists, preprocessing options,
%   and pipeline flags are defined here.
%
% OUTPUT:
%   cfg  - Configuration structure (see sections below for all fields)
%
% USAGE:
%   cfg = cts_config();
%
% =========================================================================

%% 1. ANALYSIS SELECTION
% -------------------------------------------------------------------------
% Enable one or more analysis types. Each enabled analysis will run in
% sequence for every subject/trial/condition.

cfg.analysis.coherence  = true;     % Cortical coherence (CM toolbox)
cfg.analysis.TRF        = false;    % Temporal Response Function (mTRF toolbox)
cfg.analysis.ERP        = false;    % Event-Related Potential / Evoked response

%% 2. ANALYSIS SPACE
% -------------------------------------------------------------------------
cfg.space.sensor = false;
cfg.space.source = true;

%% 3. PERMUTATION STATISTICS
% -------------------------------------------------------------------------
cfg.perm_stat = false;              % true = run surrogate/permutation pass as well

%% 4. SUBJECTS
% -------------------------------------------------------------------------
cfg.subjects.exclude_ica = {        % ICA not corrected - will be skipped
    'Meg3925',...
    'Meg3949'
};

cfg.subjects.exclude_diagnosis = {  % Diagnostic / IQ concerns - will be skipped
    'Meg6666'
};

%% 5. REFERENCE SIGNALS
% -------------------------------------------------------------------------
% Define which signals to use as reference in the analysis.
% Each row: {label, type, file_key}
%
% type options:
%   'audio'     -> envelope extracted automatically via
%   envelope_extraction()
%   'lips'      -> lip aperture signal from video landmarks (lips_apperture())
%   'envelope'  -> signal already in low-frequency format, used as-is
% 
% file_key options (for 'audio' type):
%   'snd_global'    -> global (mixed) audio
%   'snd_att'       -> attended speech
%   'snd_noise'     -> noise stream only

cfg.ref_sources = {
%    label                type       file_key
    'global sound',      'audio',   'snd_global';
    'attended speech',   'audio',   'snd_att';
    'noise',             'audio',   'snd_noise';    % uncomment to include
    'mouth aperture',    'lips',    ''
};

%% 6. PREPROCESSING OPTIONS
% -------------------------------------------------------------------------

% --- Bandpass filter (applied to MEG/EEG signal before analysis) ---
cfg.preproc.bandpass.do         = false;
cfg.preproc.bandpass.freq_low   = 1;     % High-pass cutoff (Hz)
cfg.preproc.bandpass.freq_high  = 40;    % Low-pass cutoff (Hz)
cfg.preproc.bandpass.trans_low  = 0.5;   % Transition bandwidth low (Hz)
cfg.preproc.bandpass.trans_high = 5;     % Transition bandwidth high (Hz)

% --- EEG rereferencing (only applied to .cnt files) ---
cfg.preproc.reref.do            = false;
cfg.preproc.reref.type          = 'average';  % 'average' | 'linked_mastoids' | 'none'

% --- Artifact rejection (amplitude threshold) ---
cfg.preproc.artifact.do         = false;
cfg.preproc.artifact.threshold  = 3000e-15;  % Amplitude threshold (T for MEG, V for EEG)
cfg.preproc.artifact.margin_sec = 0.2;       % Margin around detected artifacts (seconds)

% --- Z-score normalisation (applied per channel before analysis) ---
cfg.preproc.zscore.do           = false;
cfg.preproc.zscore.baseline_sec = [];    % [] = use full signal; [t1 t2] = baseline window

%% 7. COHERENCE PARAMETERS
% -------------------------------------------------------------------------
cfg.cm.quantum      = 2;        % Analysis window half-length (seconds)
cfg.cm.window_type  = 'boxcar'; % Spectral window ('boxcar', 'hanning', 'hamming')
cfg.cm.freq_max     = 20;       % Maximum frequency of interest for CM.Find (Hz)

% --- Epoch count equalisation across conditions ---
% Set equalize_nave to true if trial durations differ between conditions
% and you want to match the number of epochs used in each.
% equalize_vids: M x 2 matrix [vid_index, threshold_cond]
%   threshold > 0 : equalise when n_cond <= threshold
%   threshold < 0 : equalise when n_cond > abs(threshold)
% n_rm: vector (1 x Ncond) of epochs to remove per condition
cfg.cm.equalize_nave  = false;
cfg.cm.n_rm           = [];      % e.g. round(mean(nave(:,:,1:47),3)-330)'
cfg.cm.equalize_vids  = [2,  4;  % vid2: equalise if n_cond <= 4
                         3, -4]; % vid3: equalise if n_cond > 4

%% 8. TRF PRAMETERS
% -------------------------------------------------------------------------
cfg.trf.n_folds       = 10;
cfg.trf.refs_keep     = [1 2 3 4];
cfg.trf.att_ref_index = 2;

cfg.trf.cond_sets = {[1 2] [3 4] [5 6] [7 8] [9 10] ...
                     [1 3] [2 4] [5 7] [6 8] [9 10]};

% Forward TRF bands — columns: label | band | width | ds | t_pre | t_post | t_buff | lambda
trf_fwd = {
%  label               band        width     ds   t_pre  t_post  t_buff  lambda
  'phrasal_broad',   [0.2 1.5],  [0.2 0.5],  50,  -1000,  2000,   500,   2.^(-10:0);
  'syllabic',        [2   8  ],  [1   2  ],  10,   -100,   500,   100,   2.^(-10:0);
  'phrasal_narrow',  [0.1 0.5],  [0.1 0.2],  50,  -1000,  2000,   500,   2.^(-10:0);
  'syllabic_narrow', [2   4  ],  [0.5 1  ],  20,   -500,  1000,   200,   2.^(-10:0);
  'delta_theta',     [0.2 4  ],  [0.1 1  ],  20,  -1000,  2000,   500,   2.^(-10:0);
  'theta',           [4   8  ],  [1   2  ],  10,   -100,   500,   100,   2.^(-10:0);
};

fields = {'label','band','width','ds','t_pre','t_post','t_buff','lambda'};
for i = 1:size(trf_fwd, 1)
    for j = 1:length(fields)
        cfg.trf.fwd_bands(i).(fields{j}) = trf_fwd{i,j};
    end
end


%% 10. BEAMFORMING / SOURCE SPACE PARAMETERS
% -------------------------------------------------------------------------
cfg.beam.n_sensors    = 306;    % Number of MEG sensors
cfg.beam.ref_index    = 2;      % Reference signal index for CSD (2 = attended speech)
cfg.beam.dxyz         = 1;      % Source grid spacing (mm)
cfg.beam.fold_fs      = 'MSI-CTS'; % FreeSurfer folder name for overlays
 
% Forward solution path — use {sub_name} as placeholder for subject name
cfg.beam.fwd_file     = '/path/to/subjects/{sub_name}/meg/{sub_name}_from_MNI-MSI-CTS-5-src-fwd.fif';
cfg.beam.subjects_dir = '/path/to/freesurfer/subjects';
cfg.beam.group_fold   = '/path/to/subjects/group_CTS';
 
% Frequency bands for source projection
% Each entry needs: .label (string) and .Find (frequency indices into CM.f)
cfg.beam.freq_bands(1).label = 'phrasal';   % ~1 Hz
cfg.beam.freq_bands(1).Find  = 2;
cfg.beam.freq_bands(2).label = 'syllabic';  % ~4-8 Hz
cfg.beam.freq_bands(2).Find  = 9:17;

%% 11. FREQUENCY PARAMETERS (envelope extraction)
% -------------------------------------------------------------------------
cfg.freq.mel_low    = 150;      % MEL filterbank lower bound (Hz)
cfg.freq.mel_high   = 7000;     % MEL filterbank upper bound (Hz)
cfg.freq.mel_bands  = 31;       % Number of MEL bands

cfg.freq.align_bp   = [50 330]; % Bandpass range used for alignment verification (Hz)

%% . CONDITIONS
% -------------------------------------------------------------------------
% Each condition is defined by a label and a target vector encoding
% the experimental state: [visual_type; noise_block1; noise_block2; noise_block3]
%
% visual_type: 1 = video  |  2 = audio-only  |  3 = picture
% noise flags: 0 = clean  |  1 = noise

cfg.conditions(1).label  = 'video, clean';
cfg.conditions(1).target = [0; 0; 0; 1];

cfg.conditions(2).label  = 'video, noise';
cfg.conditions(2).target = [1; 1; 1; 1];

cfg.conditions(3).label  = 'picture, clean';
cfg.conditions(3).target = [3; 0; 0; 0];

cfg.conditions(4).label  = 'picture, noise';
cfg.conditions(4).target = [3; 1; 1; 1];

cfg.conditions(5).label  = 'audio-only, clean';
cfg.conditions(5).target = [2; 0; 0; 0];

cfg.conditions(6).label  = 'audio-only, noise';
cfg.conditions(6).target = [2; 1; 1; 1];

%% . PIPELINE BEHAVIOUR
% -------------------------------------------------------------------------
cfg.expected_trials   = 4;      % Expected number of trials per subject
cfg.overwrite_results = false;  % false = skip subject if result file already exists
cfg.verbose           = true;   % true = print progress to console
