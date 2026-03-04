function MISC = load_MISC(file, perm_stat)


MISC = struct('signal', [], 'Fs', []);
[~, ~, extension] = fileparts(file);

if strcmp(extension, '.fif')
    raw = fiff_setup_read_raw(file);
    MISC.Fs = double(raw.info.sfreq);
    MISC.sensors = get_sensors(raw);
    n_samples = raw.last_samp - raw.first_samp + 1;

    if perm_stat
        MISC.signal = zeros(1, n_samples);
    else
        picks_misc = fiff_pick_channels(raw.info.ch_names,{'MISC007'});
        MISC.signal = fiff_read_raw_segment(raw, 0, inf, picks_misc);
    end

elseif strcmp(extension, '.cnt')
    cfg = [];
    cfg.dataset = file;
    data = ft_preprocessing(cfg);
    MISC.Fs = double(data.fsample);
    picks_misc = find(strcmp(data.label, {'s2_E258'}));

    if perm_stat
        MISC.signal = zeros(1, size(data.trial{1}, 2));
    else

        MISC.signal = data.trial{1}(picks_misc, :);
    end
else
    error('File type not supported: should be .fif or .cnt')
end

MISC.signal = MISC.signal(:)';
end