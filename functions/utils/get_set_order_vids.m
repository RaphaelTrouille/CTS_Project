function [set, order, vids] = get_set_order_vids(subfold, expected_n_vids)
% GET_SET_ORDER_VIDS  Extract stimulus set, order, and video IDs from a
%                     subject's M/EEG .fif or .cnt filenames using regex-based parsing.
%
% DESCRIPTION:
%   Scans all .fif/.cnt files in a subject folder and parses their names to
%   extract the stimulus set identifier, presentation order, and video IDs.
%   Parsing is regex-based to handle flexible and evolving naming conventions.
%   If no files match the primary pattern (set + order + vid), a fallback
%   pattern (vid only) is attempted, with set inferred from fixed groupings.
%   Final outputs are validated against expected counts.
%
% INPUTS:
%   subfold         - Path to the subject's MEG/EEG data folder (.fif or .cnt files)
%   expected_n_vids - (Optional) Expected number of unique video IDs.
%                     Default: 4. Set to 0 to disable this validation check.
%
% OUTPUTS:
%   set    - Unique stimulus set identifier(s) found — numeric or string
%            (e.g. 1 or 'A'). Should resolve to a scalar after validation.
%   order  - Unique presentation order number(s). Should resolve to a scalar.
%   vids   - Vector of unique video IDs (integers) for this subject.
%
% SUPPORTED NAMING PATTERNS (regex-based, order-insensitive):
%   Primary   : requires 'set', 'order'/'ord', and 'vid' tokens
%     e.g. 'meg01_set1_order2_vid09_tsss_mc_ica.fif'
%     e.g. 'meg01_setA_ord2_vid3.fif'
%     e.g. 'sub-01_set2_ord1_vid12.fif'
%   Fallback  : vid only (set inferred from groupings defined in the code)
%     e.g. 'meg01_vid12_tsss_mc_ica.fif'
%
% CUSTOMIZATION:
%   - To support new naming formats, update the regex patterns in the
%     "REGEX PATTERNS" section below.
%   - To change set inference groupings in the fallback, update the
%     SET_GROUPS cell array in the "FALLBACK SET INFERENCE" section.
%
% NOTES:
%   - Set identifiers can be numeric (e.g. 1, 2) or alphabetic (e.g. A, B, C).
%   - Video numbers can be 1 or 2 digits, zero-padded or not.
%   - Errors are raised with informative messages if validation fails.
%
% USAGE:
%   [set, order, vids] = get_set_order_vids(sub_fold);
%   [set, order, vids] = get_set_order_vids(sub_fold, 6);  % expect 6 videos
%   [set, order, vids] = get_set_order_vids(sub_fold, 0);  % skip vid count check
%
% -------------------------------------------------------------------------

    % --- Default arguments ---
    if nargin < 2 || isempty(expected_n_vids)
        expected_n_vids = 4;
    end

    % --- Scan .fif and .cnt files ---
    files = [dir(fullfile(subfold, '*.fif')); ...
             dir(fullfile(subfold, '*.cnt'))];
    if isempty(files)
        error('get_set_order_vids: No .fif or .cnt files found in %s', subfold);
    end
    Nfiles = length(files);

    % ---------------------------------------------------------------------
    % REGEX PATTERNS - update here to support new naming formats
    % set   : matches 'set' followed by digits or letters (e.g. set1, setA)
    % order : matches 'order' or 'ord' followed by digits (e.g. order2, ord1)
    % vid   : matches 'vid' followed by 1-2 digits        (e.g. vid3, vid09)
    % ---------------------------------------------------------------------

    pat_set   = 'set([A-Za-z0-9]+?)(?:_|\.| )';
    pat_order = '(?:order|ord)(\d+?)(?:_|\.| )';
    pat_vid   = 'vid(\d{1,2})(?:_|\.| |$)';

    % --- Preallocate with NaN / empty ---
    raw_set   = cell(1, Nfiles);
    raw_order = nan(1, Nfiles);
    raw_vids  = nan(1, Nfiles);

    % --- Primary parsing: set + order + vid ---
    for n_file = 1:Nfiles
        fname = files(n_file).name;

        tok_set   = regexp(fname, pat_set, 'tokens', 'once');
        tok_order = regexp(fname, pat_order, 'tokens', 'once');
        tok_vid   = regexp(fname, pat_vid, 'tokens', 'once');

        if ~isempty(tok_set), raw_set{n_file} = tok_set{1}; end
        if ~isempty(tok_order), raw_order(n_file) = str2double(tok_order{1}); end
        if ~isempty(tok_vid), raw_vids(n_file) = str2double(tok_vid{1}); end
    end

    % Convert set tokens to numeric if all are numeric strings
    set_vals = raw_set(~cellfun(@isempty, raw_set));
    if all(cellfun(@(s) ~isnan(str2double(s)), set_vals))
        set_vals = cellfun(@str2double, set_vals);
    end

    order = unique(raw_order(~isnan(raw_order)));
    vids  = unique(raw_vids(~isnan(raw_vids)));

    if iscell(set_vals)
        set = unique(set_vals);
    else
        set = unique(set_vals);
    end

    % --- Fallback: vid-only filenames, set inferred from groupings ---
    if isempty(vids)
        fprintf('get_set_order_vids: Primary pattern matched nothing. Trying fallback...\n');

        % -----------------------------------------------------------------
        % FALLBACK SET INFERENCE — update groupings here if sets change
        % Each row: {set_label, [vid IDs belonging to this set]}
        % -------------------------------------------------------------------
        SET_GROUPS = {
            1, [1  4  9  12];
            2, [2  6  7  11];
            3, [3  5  8  10];
        };

        raw_set   = cell(1, Nfiles);
        raw_order = nan(1, Nfiles);
        raw_vids  = nan(1, Nfiles);

        for n_file = 1:Nfiles
            fname   = files(n_file).name;
            tok_vid = regexp(fname, pat_vid, 'tokens', 'once');

            if ~isempty(tok_vid)
                v = str2double(tok_vid{1});
                raw_vids(n_file)  = v;
                raw_order(n_file) = 1;  % Default order when not in filename

                % Infer set from groupings
                for g = 1:size(SET_GROUPS, 1)
                    if any(v == SET_GROUPS{g, 2})
                        raw_set{n_file} = num2str(SET_GROUPS{g, 1});
                        break
                    end
                end
            end
        end

        set_vals = raw_set(~cellfun(@isempty, raw_set));
        if all(cellfun(@(s) ~isnan(str2double(s)), set_vals))
            set_vals = cellfun(@str2double, set_vals);
        end
        if iscell(set_vals), set = unique(set_vals);
        else,                set = unique(set_vals);
        end
        order = unique(raw_order(~isnan(raw_order)));
        vids  = unique(raw_vids(~isnan(raw_vids)));
    end

    % --- Validation ---
    if length(set) ~= 1
        error('get_set_order_vids: Expected 1 unique set, found %d: [%s]\n(subfold: %s)', ...
              length(set), strjoin(arrayfun(@num2str, set(:)', 'UniformOutput', false), ', '), subfold);
    end
    if length(order) ~= 1
        error('get_set_order_vids: Expected 1 unique order, found %d: [%s]\n(subfold: %s)', ...
              length(order), num2str(order), subfold);
    end
    if expected_n_vids > 0 && length(vids) ~= expected_n_vids
        error('get_set_order_vids: Expected %d unique video IDs, found %d: [%s]\n(subfold: %s)', ...
              expected_n_vids, length(vids), num2str(vids), subfold);
    end
end