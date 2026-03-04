function files = get_subject_files(sub_fold, snd_dir, sub_name, set_order, n_vid)
    id_vid = [set_order '_vid' num2str(n_vid)];
    files.id = id_vid;

    % Flexible localization of the final folder
    search_pattern = fullfile(snd_dir, '**', [id_vid '_att.wav']);
    found_file = dir(search_pattern);

    if ~isempty(found_file)
        actual_path = found_file(1).folder;
        files.snd_att = fullfile(actual_path, [id_vid '_att.wav']);
        files.snd_global = fullfile(actual_path, [id_vid '_global.wav']);
        files.snd_noise = fullfile(actual_path, [id_vid '_noise.wav']);
    else
        warning('Files not found for %s in %s', id_vid, snd_dir)
    end
    files.meg_file = fullfile(sub_fold, [sub_name '_' id_vid '_tsss_mc_ica.fif']);
    files.syncfile = fullfile(sub_fold, [sub_name '_' id_vid '_meg_sound_sync.mat']);
    
end