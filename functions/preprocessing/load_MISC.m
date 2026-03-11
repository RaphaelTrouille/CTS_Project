function MISC = load_MISC(file, perm_stat)
% LOAD_MISC  Load a MISC (miscellaneous) channel from a MEG recording file.
%
% DESCRIPTION:
%   Reads a single MISC channel from a M/EEG data file (.fif or .cnt format)
%   and returns it as a standardized structure containing the signal and
%   its sampling rate. In permutation testing mode, the signal is replaced
%   by a zero vector of the correct length to allow label shuffling without
%   reloading data.
%
% INPUTS:
%   file       - Full path to the MEG data file (.fif or .cnt)
%   perm_stat  - Boolean flag; if true, signal is set to zeros (permutation
%                testing mode — skips actual data reading)
%
% OUTPUTS:
%   MISC       - Structure with fields:
%                  .signal  : MISC channel signal as a row vector (1 x N)
%                  .Fs      : sampling rate in Hz
%                  .sensors : sensor layout structure (only for .fif files,
%                             via get_sensors)
%
% SUPPORTED FORMATS:
%   .fif  - Elekta/Neuromag format, read via MNE/FieldTrip FIFF routines.
%           Extracts channel 'MISC007'.
%   .cnt  - Neuroscan format, read via FieldTrip ft_preprocessing.
%           Extracts channel 's2_E258'.
%
% DEPENDENCIES:
%   - fiff_setup_read_raw, fiff_pick_channels, fiff_read_raw_segment (MNE/FieldTrip)
%   - ft_preprocessing (FieldTrip)
%   - get_sensors.m (cartographie_motrice)
%
% USAGE:
%   MISC = load_MISC('path/to/recording.fif', false);
%
% -------------------------------------------------------------------------

MISC = struct('signal', [], 'Fs', []);
[~, ~, extension] = fileparts(file);

if strcmp(extension, '.fif')

    % Read raw .fif file metadata and extract sampling rate and sensor info
    raw = fiff_setup_read_raw(file);
    MISC.Fs = double(raw.info.sfreq);
    MISC.sensors = get_sensors(raw);
    n_samples = raw.last_samp - raw.first_samp + 1;

    if perm_stat
        % Permutation mode: return a zero signal of the correct length
        MISC.signal = zeros(1, n_samples);
    else
        % Read only the MISC007 channel
        picks_misc = fiff_pick_channels(raw.info.ch_names,{'MISC007'});
        MISC.signal = fiff_read_raw_segment(raw, 0, inf, picks_misc);
    end

elseif strcmp(extension, '.cnt')

    % Load .cnt file via FieldTrip preprocessing
    cfg = [];
    cfg.dataset = file;
    data = ft_preprocessing(cfg);
    MISC.Fs = double(data.fsample);
    picks_misc = find(strcmp(data.label, {'s2_E258'}));

    if perm_stat
        % Permutation mode: return a zero signal of the correct length
        MISC.signal = zeros(1, size(data.trial{1}, 2));
    else
        % Extract the MISC channel (s2_E258)
        MISC.signal = data.trial{1}(picks_misc, :);
    end
else
    error('File type not supported: should be .fif or .cnt');
end

% Ensure signal is a row vector
MISC.signal = MISC.signal(:)';
end