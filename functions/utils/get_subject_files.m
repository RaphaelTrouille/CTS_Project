function files = get_subject_files(sub_fold, snd_dir, sub_name, set, order, n_vid)
% GET_SUBJECT_FILES  Build file paths for all data sources associated with
%                    a given subject and video stimulus.
%
% DESCRIPTION:
%   Constructs and returns a structure containing all file paths required
%   to process one trial (subject x video). Audio files are located via a
%   recursive directory search to handle flexible folder structures. MEG
%   and synchronization file paths are built directly from the subject
%   folder and naming conventions.
%
% INPUTS:
%   sub_fold   - Path to the subject's MEG data folder
%   snd_dir    - Root path to the audio stimuli directory (searched recursively)
%   sub_name   - Subject identifier string (e.g. 'meg01')
%   set_order  - Stimulus set identifier string (e.g. 'setA')
%   n_vid      - Video number (integer)
%
% OUTPUTS:
%   files  - Structure with fields:
%              .id          : trial identifier string
%              ('<set_order>_vid<n_vid>') (e.g. 'set1_order2_vid3')
%              .snd_att     : path to the attended speech WAV file
%              .snd_global  : path to the global (mixed) audio WAV file
%              .snd_noise   : path to the noise WAV file
%              .meg_file    : path to the preprocessed MEG .fif file
%              .sync_file    : path to the MEG/audio synchronization .mat file
%
% FILE NAMING CONVENTIONS:
%   Audio  : '<id>_att.wav', '<id>_global.wav', '<id>_noise.wav'
%   MEG    : '<sub_name>_<id>_tsss_mc_ica.fif'
%   Sync   : '<sub_name>_<id>_meg_sound_sync.mat'
%
% NOTES:
%   - Audio files are searched recursively under snd_dir using a glob pattern.
%     A warning is issued if no matching files are found.
%   - Only the first match is used if multiple files are found.
%
% USAGE:
%   files = get_subject_files(sub_fold, snd_dir, 'meg01', 'set1_order2', 3);
%
% -------------------------------------------------------------------------

    % Build trial identifier from set order and video number
    set_order = ['set' num2str(set) '_order' num2str(order)];
    id_vid = [set_order '_vid' num2str(n_vid)];
    files.id = id_vid;

    % Recursively search for the attended WAV file to locate the audio
    % folder
    search_pattern = fullfile(snd_dir, '**', [id_vid '_att.wav']);
    found_file = dir(search_pattern);

    if ~isempty(found_file)
        % Use the found folder to build all audio file paths
        actual_path = found_file(1).folder;
        files.snd_att = fullfile(actual_path, [id_vid '_att.wav']);
        files.snd_global = fullfile(actual_path, [id_vid '_global.wav']);
        files.snd_noise = fullfile(actual_path, [id_vid '_noise.wav']);
    else
        warning('Files not found for %s in %s', id_vid, snd_dir)
    end

    % Build MEG and synchronization file paths from subject folder
    files.meg_file = fullfile(sub_fold, [sub_name '_' id_vid '_tsss_mc_ica.fif']);
    files.sync_file = fullfile(sub_fold, [sub_name '_' id_vid '_meg_sound_sync.mat']);
    
end