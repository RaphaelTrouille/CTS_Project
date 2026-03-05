function [CM, bad] = prepare_CM_data(CM, dec, L, file)
% PREPARE_CM_DATA  Load and prepare MEG/EEG data into a CM structure for
%                  cortical mapping analysis.
%
% DESCRIPTION:
%   Reads a raw MEG (.fif) or EEG (.cnt) data file, defines the analysis
%   time window based on a temporal offset (dec) and length (L), detects
%   artefacts (flat channels, boundary periods), and populates the CM
%   structure with channel labels, picks, and metadata.
%
% INPUTS:
%   CM    - CM structure to populate (must already exist)
%   dec   - Temporal offset in samples (from realignment). Positive = MEG leads audio.
%   L     - Desired analysis segment length in samples
%   file  - Full path to the raw data file (.fif or .cnt)
%
% OUTPUTS:
%   CM    - Updated CM structure with fields:
%             .picksMEEG : channel indices used for analysis
%             .label     : channel name list
%             .infile    : path to the source file
%             .CSD       : (fif only) initialized to [] for CSD computation
%   bad   - Binary artefact mask (1 x N samples); 1 = artefact, 0 = clean.
%           Includes: pre/post-analysis margins, flat channels, and a 2s
%           smoothing buffer around all flagged regions.
%
% SUPPORTED FORMATS:
%   .cnt  - Neuroscan EEG, read via FieldTrip ft_preprocessing.
%           Uses channels 1:256. Flat channels are auto-detected and flagged.
%           Signal is mean-centred across channels before storage.
%   .fif  - Elekta/Neuromag MEG, read via MNE/FieldTrip FIFF routines.
%           Bad channel mask is read from a dedicated 'bad' channel.
%
% DEPENDENCIES:
%   - ft_preprocessing        (FieldTrip)
%   - fiff_setup_read_raw     (MNE/FieldTrip)
%   - fiff_read_raw_segment   (MNE/FieldTrip)
%   - fiff_pick_channels      (MNE/FieldTrip)
%   - get_sensors.m           (custom, must be on MATLAB path)
%
% USAGE:
%   [CM, bad] = prepare_CM_data(CM, dec, L, 'path/to/data.fif');
%
% -------------------------------------------------------------------------

    [~, ~, extension] = fileparts(file);
    
    % Format detection and data loading
    if strcmp(extension, '.cnt')
        % --- EEG: FieldTrip preprocessing ---
        cfg         = [];
        cfg.dataset = file;
        data = ft_preprocessing(cfg);
        Fs = double(data.fsample);
        % Define analysis window with temporal offset
        tdeb = data.sampleinfo(1) + max(dec, 0);
        tfin = min(tdeb + L - round(2*Fs), data.sampleinfo(2));
    
        % --- Artefact mask: flag flat pre/post-analysis margins ---
        bad = zeros(1, length(data.trial{1}));
        bad(1:min(end, tdeb)) = 1;
        bad(max(1, tfin):end) = 1;
    
        % --- Artefact mask: flag flat (dead) channels ---
        bad(max(data.trial{1}(1:256,:),[],1) == 0) = 1;
    
        % --- Smooth artefact mask with a 2s buffer on each side ---
        bad = double(conv(double(bad), ones(1, round(2*Fs)), 'same') > 0.5);
    
        % --- Mean-centre each sample across the 256 EEG channels ---
        data.trial{1}(1:256,:) = bsxfun(@minus, data.trial{1}(1:256,:), ...
            mean(data.trial{1}(1:256,:), 1));
    
        % --- Populate CM structure ---
        CM.picksMEEG = 1:256;
        CM.label = data.label;
        CM.infile = file;
    
    elseif strcmp(extension, '.fif')
        % --- MEG: MNE/FieldTrip FIFF routines ---
        raw = fiff_setup_read_raw(file);
        Fs = double(raw.info.sfreq);
    
        % Define analysis window with temporal offset
        tdeb = double(raw.first_samp) + max(dec, 0);
        tfin = min(tdeb+L, double(raw.last_samp));
    
        % Read bad channel mask from dedicated 'bad' channel
        bad = fiff_read_raw_segment(raw, 0, inf, ...
            fiff_pick_channels(raw.info.ch_names,{'bad'}));
    
        % --- Populate CM structure ---
        sensors     = get_sensors(raw);
        CM.CSD      = [];   % Placeholder for Current Source Density
        CM.label    = raw.info.ch_names(sensors.picksMEG);
        CM.infile   = file;
    
    else
        error('prepare_CM_data: Unsupported file format "%s". Expected .fif or .cnt.', extension)
    end
end