function [meg_dir, subjects, deriv_dir, snd_dir, vid_dir] = setup_environment()
% SETUP_ENVIRONMENT Initialize the project environment for the CTS
% experiment.
%
% DESCRIPTION:
%   This function sets up all necessary paths and directories for the
%   SpeechTrack MEG experiment pipeline. It adds project functions, main
%   scripts, and required toolboxes (FieldTrip, SM8, AMToolbox, etc.) to
%   the MATLAB path. It also definesdata-related directories and retrieves
%   the list of subjects from the MEG data folder.
%
% OUTPUTS: 
%   meg_dir    - Path to the raw MEG data directory
%   subjects   - Structure array of subject folders (from dir())
%   deriv_dir  - Path to the derivatives/processed data directory
%   snd_dir    - Path to the audio stimuli directory
%   vid_dir    - Path to the video stimuli directory
%
% DEPENDENCIES:
%   - FieldTrip (v. 20161113)
%   - SPM8
%   - AMToolbox
%   - cartographie_motrice
%   - ARfit
%   - LocComp
%   - gap_statistics
%
% USAGE:
%   [meg_dir, subjects, deriv_dir, snd_dir, vid_dir] = setup_environment();
%
% NOTES:
%   - Toolbox paths are defined as absolute paths and must be adapted
%     to your local environment (tb_dir variable).
%   - MEG data is expected on an external drive (/Volumes/Elements/...).
%     A warning is issued if the folder is not found.
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

    % 1. Root of the project based script path
    project_root = fileparts(mfilename('fullpath'));
    
    % 2. Add functions/methods/scripts relative to the path
    addpath(genpath(fullfile(project_root, 'functions')));
    addpath(fullfile(project_root, 'main_scripts'));
    
    % 3. Toolboxes
    % Absolute path to change according to your environment
    tb_dir = '/Users/raphaeltrouille/Desktop/toolboxes';
    to_add = {'cartographie_motrice'...
    'cartographie_motrice/CM_extra'...
    'spm8'...
    'arfit'...
    'fieldtrip-20161113'...
    'fieldtrip-20161113/utilities'...
    'fieldtrip-20161113/src'...
    'fieldtrip-20161113/external/mne' ...
    'fieldtrip-20161113/external/fastica'...
    'fieldtrip-20161113/external/egi_mff/'...
    'fieldtrip-20161113/fileio'...
    'fieldtrip-20161113/forward'...
    'fieldtrip-20161113/plotting'...
    'amtoolbox/code/general' ...
    'amtoolbox/code/thirdparty/ltfat'...
    'amtoolbox/code/thirdparty/ltfat/auditory'...
    'LocComp'...
    'gap_statistics'};
    
    for n_add = 1:length(to_add)
        addpath(fullfile(tb_dir,to_add{n_add}));
    end
    %ft_defaults; % Very important for FieldTrip
    
    
    % 4. Data related paths
    
    %Internal
    %meg_dir = fullfile(project_root, 'data', 'raw', 'meg');
    %snd_dir = fullfile(project_root, 'data', 'audio', 'meg');
    %vid_dir = fullfile(project_root, 'data', 'video');
    
    %External
    external_root = '/Volumes/Elements/expe_SpeechTrack';
    meg_dir = fullfile(external_root, 'meg');
    snd_dir = fullfile(external_root, 'ready_stim');
    vid_dir = fullfile(external_root, 'clean_vids');
    deriv_dir = fullfile(project_root, 'data', 'raw', 'meg');
    
    if ~exist(meg_dir, 'dir')
        warning('MEG data folder not found: %s', meg_dir);
        subjects = [];
    else
        subjects = dir(fullfile(meg_dir, 'meg*'));
        fprintf('>>> %d subjects found.\n', length(subjects));
    end
    
    disp('>>> CTS environment ready. Functions and Toolboxes loaded.');
end