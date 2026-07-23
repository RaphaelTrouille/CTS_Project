%% build_group_cohgrad.m
% Builds ONCE two large consolidated .mat files (real data + permutation)
% from individual subject files to avoid reloading 144+ files for every analysis.
%
% Output:
%   group_cohgrad_all_subjects.mat       (real data)
%   group_cohgrad_all_subjects_perm.mat  (permuted data / null distribution)
%
% Dimensions of all_cohgrad: (subjects x pairs x frequencies x conditions x references)
%   e.g., 144 x 102 x 21 x 9 x 4

%% Init
cfg = cts_config;
[meg_dir, ~, deriv_dir, ~, ~] = setup_environment;

matfiles_dir = '/mnt/usb-HardDrive_MB/expe_SpeechTrack/matfiles';

n_subjects = length(cfg.subjects.include);
n_pairs    = 102;
n_freq     = 21;   % Adjust if needed, verified via freq_axis below
n_cond     = 10;    % Adjust based on the actual number of conditions in CMall
n_ref      = 4;

%% Build both versions (real + perm) using the same function
fprintf('=== Building REAL data ===\n');
real_suffix = '_main_coh.mat';
[all_cohgrad, subject_ids, subject_valid, subject_nave, ...
 freq_axis, cm_labels, g1, g2] = load_all_subjects( ...
    cfg, matfiles_dir, real_suffix, n_subjects, n_pairs, n_freq, n_cond, n_ref);

fprintf('\n=== Building PERMUTED data (null distribution) ===\n');
perm_suffix = '_main_coh_perm_stat.mat';
[all_cohgrad_perm, subject_ids_perm, subject_valid_perm, subject_nave_perm, ...
 freq_axis_perm, cm_labels_perm, g1_perm, g2_perm] = load_all_subjects( ...
    cfg, matfiles_dir, perm_suffix, n_subjects, n_pairs, n_freq, n_cond, n_ref);

%% Sanity check: matching metadata between real and perm (should always match)
if ~isequal(freq_axis, freq_axis_perm)
    warning('freq_axis differs between real and perm data!')
end
if ~isequal(cm_labels, cm_labels_perm)
    warning('cm_labels differs between real and perm data!')
end
if ~isequal(g1, g1_perm) || ~isequal(g2, g2_perm)
    warning('g1/g2 differs between real and perm data!')
end

%% Log missing subjects
missing_real = subject_ids(~subject_valid);
missing_perm = subject_ids_perm(~subject_valid_perm);

log_dir = fullfile(matfiles_dir, 'logs');
if ~exist(log_dir, "dir"), mkdir(log_dir); end

if ~isempty(missing_real)
    writecell(missing_real, fullfile(log_dir, 'missing_subjects_real.txt'));
    fprintf('%d missing subjects (real), see missing_subjects_real.txt\n', length(missing_real));
end
if ~isempty(missing_perm)
    writecell(missing_perm, fullfile(log_dir, 'missing_subjects_perm.txt'));
    fprintf('%d missing subjects (perm), see missing_subjects_perm.txt\n', length(missing_perm));
end

%% Save - real data
out_real = fullfile(matfiles_dir, 'group_cohgrad_all_subjects.mat');
save(out_real, ...
    'all_cohgrad', 'subject_ids', 'subject_valid', 'subject_nave', ...
    'freq_axis', 'cm_labels', 'g1', 'g2', '-v7.3');
fprintf('\nSaved: %s\n', out_real);

%% Save - permuted data
out_perm = fullfile(matfiles_dir, 'group_cohgrad_all_subjects_perm.mat');
save(out_perm, ...
    'all_cohgrad_perm' , 'subject_ids_perm', 'subject_valid_perm', 'subject_nave_perm', ...
    'freq_axis_perm', 'cm_labels_perm', 'g1_perm', 'g2_perm', '-v7.3');
fprintf('Saved: %s\n', out_perm);

fprintf('\nFinished. %d/%d valid subjects (real), %d/%d valid subjects (perm).\n', ...
    sum(subject_valid), n_subjects, sum(subject_valid_perm), n_subjects);


%% ======================================================================
function [all_cohgrad, subject_ids, subject_valid, subject_nave, ...
          freq_axis, cm_labels, g1, g2] = load_all_subjects( ...
          cfg, matfiles_dir, suffix, n_subjects, n_pairs, n_freq, n_cond, n_ref)
% Loads all subjects for a given suffix (real or perm) and stacks
% the cohgrad matrices into a single 5D array.

    all_cohgrad   = nan(n_subjects, n_pairs, n_freq, n_cond, n_ref);
    subject_ids   = cell(n_subjects, 1);
    subject_nave  = nan(n_subjects, n_cond);
    subject_valid = false(n_subjects, 1);

    freq_axis = [];
    cm_labels = {};
    g1 = [];
    g2 = [];
    first_loaded = false;

    for i = 1:n_subjects
        subject = cfg.subjects.include{i};
        subject_ids{i} = subject;

        matfile = fullfile(matfiles_dir, subject, [subject suffix]);

        if ~exist(matfile, "file")
            fprintf('  [%3d/%3d] %s: file not found (%s)\n', i, n_subjects, subject, suffix);
            continue
        end

        try
            S = load(matfile, 'CMall');
            CMall = S.CMall;
        catch ME
            fprintf('  [%3d/%3d] %s: error loading file (%s)\n', i, n_subjects, subject, ME.message);
            continue
        end

        if numel(CMall) < n_cond
            fprintf('  [%3d/%3d] %s: only %d conditions found (expected %d), subject skipped\n', ...
                i, n_subjects, subject, numel(CMall), n_cond);
            continue
        end

        ok = true;
        for c = 1:n_cond
            cg = CMall(c).cohgrad;   % Expected: (n_pairs x n_freq x n_ref)
            if ~isequal(size(cg), [n_pairs, n_freq, n_ref])
                fprintf('  [%3d/%3d] %s: cohgrad cond %d has unexpected dimensions [%s], subject skipped\n', ...
                    i, n_subjects, subject, c, num2str(size(cg)));
                ok = false;
                break
            end
            all_cohgrad(i, :, :, c, :) = cg;
            subject_nave(i, c) = double(CMall(c).nave);
        end

        if ~ok
            all_cohgrad(i, :, :, :, :) = nan;   % Reset partial data if check failed
            continue
        end

        subject_valid(i) = true;

        % Shared metadata: extracted from the first valid subject, then 
        % cross-checked for subsequent subjects
        if ~first_loaded
            freq_axis = CMall(1).f;
            cm_labels = CMall(1).label;
            g1        = double(CMall(1).g1);
            g2        = double(CMall(1).g2);
            first_loaded = true;
        else
            if ~isequal(CMall(1).f, freq_axis)
                warning('Subject %s: freq_axis differs from the reference!', subject);
            end
            if ~isequal(CMall(1).label, cm_labels)
                warning('Subject %s: channel labels differ from the reference!', subject);
            end
        end

        fprintf('  [%3d/%3d] %s: OK\n', i, n_subjects, subject);
    end
end