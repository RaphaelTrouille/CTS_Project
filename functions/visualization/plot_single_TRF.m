function fig = plot_single_TRF(CMfwd, cfg, subID)
% PLOT_SINGLE_TRF  Visualise TRF kernels and sensor-level prediction performance.
%
%   FIG = PLOT_SINGLE_TRF(CMfwd, CFG, SUBID) produces a two-panel figure for
%   one (condition × frequency band × hemisphere) entry of the CMfwd result
%   array:
%       Panel 1 – TRF weight kernels (all sensors in grey, mean in red).
%       Panel 2 – Topographic map of mean cross-validated prediction
%                 performance (r), with a bar-plot fallback when the topoplot
%                 function is unavailable or raises an error.
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
%   CMfwd   Struct array [nConds × nBands × nSides].  Each element must
%           contain at minimum the fields:
%               .rval   [nSensors × nRefs × nFolds]  Pearson r per fold.
%               .w      {nRefs}  TRF weight cell array; each cell holds a
%                       matrix [nSensors × nFolds].
%               .t      [1 × nLags]  Time-lag axis (ms).
%
%   cfg     Scalar struct with the following fields:
%
%       REQUIRED
%           .cond   (integer) Index into the first dimension of CMfwd.
%           .band   (integer) Index into the second dimension of CMfwd.
%           .side   (integer) Index into the third dimension of CMfwd.
%           .ref    (integer) Reference (envelope) index for rval and w.
%
%       OPTIONAL DISPLAY LABELS  (cell arrays of char/string)
%           .cond_labels    Labels for each condition.
%           .band_labels    Labels for each frequency band.
%           .side_labels    Labels for each hemisphere side.
%           .ref_labels     Labels for each reference signal.
%
%       OPTIONAL BEHAVIOUR
%           .save           (logical, default false)  Save the figure to disk.
%           .save_dir       (char, default './results/figures')  Output folder.
%           .show_console   (logical, default true)   Print summary to console.
%           .topoplot_fun   (function handle, default @eeg_topoplot)
%                           Topoplot function with signature f(values, clim).
%
%   subID   Character vector or string.  Subject identifier used in the
%           figure title and (optionally) the saved filename.
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
%   fig     Handle to the generated figure.
%
% -------------------------------------------------------------------------
% EXAMPLE
% -------------------------------------------------------------------------
%   cfg.cond = 1;  cfg.band = 2;  cfg.side = 1;  cfg.ref = 1;
%   cfg.cond_labels = {'Quiet','Noise'};
%   cfg.band_labels = {'Delta','Theta'};
%   cfg.side_labels = {'Left','Right'};
%   cfg.ref_labels  = {'Attended','Unattended'};
%   fig = plot_single_TRF(CMfwd, cfg, 'sub-01');
%
% =========================================================================

% =========================================================================
% 0 – Input validation and defaults
% =========================================================================

validateattributes(CMfwd, {'struct'}, {'nonempty'}, mfilename, 'CMfwd', 1);
validateattributes(cfg,   {'struct'}, {'nonempty'}, mfilename, 'cfg',   2);

if ~(ischar(subID) || isstring(subID))
    error('%s: subID must be a character vector or string.', mfilename);
end
subID = char(subID);

% --- Verify that all required cfg fields are present --------------------
requiredFields = {'cond', 'band', 'side', 'ref'};
for iF = 1:numel(requiredFields)
    if ~isfield(cfg, requiredFields{iF})
        error('%s: Missing required cfg field: cfg.%s', mfilename, requiredFields{iF});
    end
end

validateattributes(cfg.cond, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'cfg.cond');
validateattributes(cfg.band, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'cfg.band');
validateattributes(cfg.side, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'cfg.side');
validateattributes(cfg.ref,  {'numeric'}, {'scalar','integer','positive'}, mfilename, 'cfg.ref');

% --- Optional fields: behaviour -----------------------------------------
if ~isfield(cfg, 'save')         || isempty(cfg.save),         cfg.save         = false;                              end
if ~isfield(cfg, 'save_dir')     || isempty(cfg.save_dir),     cfg.save_dir     = fullfile(pwd,'results','figures');  end
if ~isfield(cfg, 'show_console') || isempty(cfg.show_console), cfg.show_console = true;                               end
if ~isfield(cfg, 'topoplot_fun') || isempty(cfg.topoplot_fun), cfg.topoplot_fun = @eeg_topoplot;                      end

% --- Convenience aliases ------------------------------------------------
c   = cfg.cond;
b   = cfg.band;
s   = cfg.side;
ref = cfg.ref;

% --- Bounds check on CMfwd dimensions -----------------------------------
if c > size(CMfwd, 1)
    error('%s: cfg.cond (%d) exceeds CMfwd condition dimension (%d).', mfilename, c, size(CMfwd,1));
end
if b > size(CMfwd, 2)
    error('%s: cfg.band (%d) exceeds CMfwd band dimension (%d).',      mfilename, b, size(CMfwd,2));
end
if s > size(CMfwd, 3)
    error('%s: cfg.side (%d) exceeds CMfwd side dimension (%d).',      mfilename, s, size(CMfwd,3));
end

% --- Resolve display labels (fall back to "Dim N" strings) --------------
condLabel = get_cfg_label(cfg, 'cond_labels', c,   sprintf('Cond %d', c));
bandLabel = get_cfg_label(cfg, 'band_labels', b,   sprintf('Band %d', b));
sideLabel = get_cfg_label(cfg, 'side_labels', s,   sprintf('Side %d', s));
refLabel  = get_cfg_label(cfg, 'ref_labels',  ref, sprintf('Ref %d',  ref));

% --- Select the result block of interest --------------------------------
CM = CMfwd(c, b, s);

% =========================================================================
% 1 – Extract prediction performance r
% =========================================================================

if ~isfield(CM, 'rval') || isempty(CM.rval)
    error('%s: CMfwd(%d,%d,%d).rval is missing or empty.', mfilename, c, b, s);
end

% rval : [nSensors × nRefs × nFolds]
% r_mean : [nSensors × 1]  – average r across cross-validation folds
rval   = CM.rval;
r_mean = extract_rval(rval, ref);

% =========================================================================
% 2 – Extract TRF weight kernels
% =========================================================================

if ~isfield(CM, 'w') || isempty(CM.w) || isempty(CM.w{ref})
    error('%s: CMfwd(%d,%d,%d).w{%d} is missing or empty.', mfilename, c, b, s, ref);
end
if ~isfield(CM, 't') || isempty(CM.t)
    error('%s: CMfwd(%d,%d,%d).t is missing or empty.', mfilename, c, b, s);
end

t      = CM.t;                      % [1 × nLags]  time-lag axis (ms)
w_raw  = squeeze(CM.w{ref});        % [nSensors × nFolds]  per-fold kernels
w_mean = mean(w_raw, 2);            % [nSensors × 1]       mean kernel

% =========================================================================
% 3 – Summary metrics (optionally printed to console)
% =========================================================================

mean_r          = mean(r_mean, 'omitnan');
max_r           = max(r_mean);
[~, bestSensor] = max(r_mean);

if cfg.show_console
    fprintf('\n[plot_single_TRF] %s | %s | %s | %s | %s\n', ...
        subID, condLabel, bandLabel, sideLabel, refLabel);
    fprintf('  Mean r (across sensors) : %.4f\n', mean_r);
    fprintf('  Max  r                  : %.4f  (sensor %d)\n', max_r, bestSensor);
end

% =========================================================================
% 4 – Figure layout
% =========================================================================

fig = figure( ...
    'Name',     ['TRF_' subID], ...
    'Color',    'w', ...
    'Position', [100 100 1250 450]);

% -------------------------------------------------------------------------
% Panel 1 – TRF weight kernels over time
% -------------------------------------------------------------------------
subplot(1, 2, 1);

% Individual sensor kernels in light grey
plot(t, w_raw,  'Color', [.85 .85 .85], 'LineWidth', 0.8);
hold on;

% Grand mean across sensors in red
plot(t, w_mean, 'r', 'LineWidth', 2.2);

% Vertical marker at zero lag
xline(0, '--k', 'LineWidth', 1);

xlabel('Time lag (ms)');
ylabel('TRF weight (a.u.)');
title(sprintf('TRF kernels  |  %s  |  %s  |  %s  |  %s', ...
    condLabel, bandLabel, sideLabel, refLabel), ...
    'Interpreter', 'none');
legend({'Sensors', 'Mean across sensors'}, 'Location', 'best');
legend boxoff;
grid on;
axis tight;
box off;

% -------------------------------------------------------------------------
% Panel 2 – Topographic map of prediction performance (r)
% -------------------------------------------------------------------------
subplot(1, 2, 2);

r_mean = r_mean(:);   % ensure column vector

% Compute colour limits; guard against all-NaN or constant-value vectors
validVals = r_mean(~isnan(r_mean));
if isempty(validVals)
    clim = [-1 1] * 1e-6;
else
    clim = [min(validVals) max(validVals)];
    if diff(clim) == 0
        clim = clim + [-1 1] * 1e-6;   % prevent flat colourmap
    end
end

try
    % Primary path: user-supplied (or default) topoplot function
    cfg.topoplot_fun(r_mean, clim);
    colorbar;
    title(sprintf('Prediction performance (r)  |  %s  |  %s', sideLabel, refLabel), ...
        'Interpreter', 'none');

catch ME
    % Fallback: simple bar plot when topoplot is unavailable
    warning('%s: Topoplot failed (%s). Using bar-plot fallback.', mfilename, ME.message);
    bar(r_mean);
    xlabel('Sensor index');
    ylabel('Prediction performance (r)');
    title(sprintf('Prediction performance (r)  |  %s  |  %s', sideLabel, refLabel), ...
        'Interpreter', 'none');
    grid on;
    box off;
end

% =========================================================================
% 5 – Optional save
% =========================================================================

if cfg.save
    if ~isfolder(cfg.save_dir)
        mkdir(cfg.save_dir);
    end
    fname = fullfile(cfg.save_dir, ...
        sprintf('TRF_%s_cond%d_band%d_side%d_ref%d.png', subID, c, b, s, ref));
    exportgraphics(fig, fname, 'Resolution', 300);
    fprintf('[plot_single_TRF] Figure saved to: %s\n', fname);
end

end   % ── end plot_single_TRF ────────────────────────────────────────────


% =========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

function label = get_cfg_label(cfg, fieldName, idx, fallback)
% GET_CFG_LABEL  Safely retrieve a display label from a cfg cell array.
%
%   LABEL = GET_CFG_LABEL(CFG, FIELDNAME, IDX, FALLBACK) returns
%   cfg.<fieldName>{IDX} if the field exists and IDX is in range;
%   otherwise returns FALLBACK.

if isfield(cfg, fieldName) && ~isempty(cfg.(fieldName))
    labels = cfg.(fieldName);
    if idx <= numel(labels)
        label = char(string(labels{idx}));
        return;
    end
end
label = fallback;

end


function r = extract_rval(rval, ref)
% EXTRACT_RVAL  Average prediction performance across cross-validation folds.
%
%   R = EXTRACT_RVAL(RVAL, REF) extracts the slice corresponding to
%   reference REF and averages over folds.
%
%   INPUT
%       rval  [nSensors × nRefs × nFolds]  per-fold Pearson r values.
%       ref   (integer) Reference index to select along dimension 2.
%
%   OUTPUT
%       r     [nSensors × 1]  Mean r across folds, NaNs omitted.

if ndims(rval) ~= 3  %#ok<ISMAT>
    error('extract_rval: expected a 3-D array [nSensors × nRefs × nFolds], got size %s.', ...
        mat2str(size(rval)));
end

nRefs = size(rval, 2);
if ref > nRefs
    error('extract_rval: cfg.ref (%d) exceeds the number of references in rval (%d).', ...
        ref, nRefs);
end

r = squeeze(rval(:, ref, :));   % [nSensors × nFolds]
r = mean(r, 2, 'omitnan');      % [nSensors × 1]

end