function [CM, bad] = prepare_CM_data(CM, dec, L, file)

[~, ~, extension] = fileparts(file);
%1. --- Type detection and data extraction ---
if strcmp(extension, '.cnt')
    % This is FieldTrip (EEG)
    cfg = [];
    cfg.dataset = file;
    data = ft_preprocessing(cfg);

    Fs = double(data.fsample);
    tdeb = data.sampleinfo(1) + max(dec, 0);
    tfin = min(tdeb + L - round(2*Fs), data.sampleinfo(2));
    % Init artefacts detection
    bad = zeros(1, length(data.trial{1}));
    bad(1:min(end, tdeb)) = 1;
    bad(max(1, tfin):end) = 1;

    % Detect flat channels
    bad(max(data.trial{1}(1:256,:),[],1) == 0) = 1;
    % Smooth 2sec over exclusion
    bad = double(conv(double(bad), ones(1, round(2*Fs)), 'same') > 0.5);
    % Centered datas
    data.trial{1}(1:256,:) = bsxfun(@minus, data.trial{1}(1:256,:), mean(data.trial{1}(1:256,:), 1));
    
    % Update CM
    CM.picksMEEG = 1:256;
    CM.label = data.label;
    CM.infile = file;

elseif strcmp(extension, '.fiff')
    raw = fiff_setup_read_raw(file);
    Fs = double(raw.info.sfreq);

    tdeb = double(raw.first_samp) + max(dec, 0);
    tfin = min(tdeb+L, double(raw.last_samp));

    CM.CSD = [];
    sensors = get_sensors(raw);
    CM.label = raw.info.ch_names(sensors.picksMEG);
    CM.infile = file;
    bad = fiff_read_raw_segment(raw, 0, inf, fiff_pick_channels(raw.info.ch_names,{'bad'}));
end