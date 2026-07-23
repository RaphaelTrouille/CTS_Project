function plot_group_TRF(rootDir, cfg_plot)

% =========================================================================
% 0) Input checks and defaults
% =========================================================================
validateattributes(rootDir, {'char', 'string'}, {'nonempty'}, mfilename, 'rootDir', 1);

validateattributes(cfg_plot, {'struct'}, {'nonempty'}, mfilename, 'cfg_plot', 2);

requiredFields = {'cond', 'band', 'ref'};
for iF = 1:numel(requiredFields)
    if ~isfield(cfg_plot, requiredFields{iF})
        error('Missing required cfg_plot field: cfg_plot.%s', requiredFields{iF});
    end
end

validateattributes(cfg_plot.cond, {'numeric'}, {'vector', 'integer', 'positive'});
validateattributes(cfg_plot.band, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(cfg_plot.side, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(cfg_plot.ref,  {'numeric'}, {'scalar', 'integer', 'positive'});

if ~isfield(cfg_plot, 'pattern') || isempty(cfg_plot.pattern)
    cfg_plot.pattern = '*.mat';
end
if ~isfield(cfg_plot, 'save') || isempty(cfg_plot.save)
    cfg_plot.save = false;
end
if ~isfield(cfg_plot, 'save_dir') || isempty(cfg_plot.save_dir)
    cfg_plot.save_dir = fullfile(rootDir, 'group_figures');
end
if ~isfield(cfg_plot, 'show_console') || isempty(cfg_plot.show_console)
    cfg_plot.show_console = false;   % bug fix: was overwriting cfg_plot.pattern
end
if ~isfield(cfg_plot, 'sensors') || isempty(cfg_plot.sensors)
    cfg_plot.sensors = [];
end
if ~isfield(cfg_plot, 'topoplot_fun') || isempty(cfg_plot.topoplot_fun)
    cfg_plot.topoplot_fun = @FT_quick_topoplot;
end

conds = cfg_plot.cond(:)';
nCond = numel(conds);   % bug fix: was numel(c)
b     = cfg_plot.band;
s     = cfg_plot.side;
ref   = cfg_plot.ref;

if ~exist(rootDir, 'dir')
    error('Input directory does not exist: %s', rootDir);
end

% =========================================================================
% 1. Find subject folders
% =========================================================================
subDirs = dir(rootDir);
subDirs = subDirs([subDirs.isdir]);
subDirs = subDirs(~ismember({subDirs.name}, {'.', '..'}));

if isempty(subDirs)
    error('No subject folders found in: %s', rootDir);
end

subjectIDs   = {};
all_r        = [];   % [nSensors x nSubjects x nCond]
all_w        = [];   % [nTime x nSensors x nSubjects x nCond]
t_ref        = [];
cfg_loaded   = [];
loadedCount  = 0;
skippedCount = 0;

% =========================================================================
% 2. Loop over subjects
% =========================================================================
for iSub = 1:numel(subDirs)

    subName = subDirs(iSub).name;
    subPath = fullfile(rootDir, subName);

    matFiles = dir(fullfile(subPath, cfg_plot.pattern));

    if isempty(matFiles)
        warning('No .mat file found in %s (pattern: %s). Skipping.', subPath, cfg_plot.pattern);
        skippedCount = skippedCount + 1;
        continue
    end

    matPath = fullfile(subPath, matFiles(1).name);

    try
        S = load(matPath);

        if ~isfield(S, 'CMall')
            warning('Missing CMall in %s. Skipping.', matFiles(1).name);
            skippedCount = skippedCount + 1;
            continue
        end
        if ~isfield(S, 'cfg')
            warning('Missing cfg in %s. Proceeding without file-level labels.', matFiles(1).name);
        end

        CMall = S.CMall;

        if isfield(S, 'cfg') && isempty(cfg_loaded)
            cfg_loaded = S.cfg;
        end

        % --- Check all conditions are present before storing anything ----
        subValid = true;
        r_sub_all = [];
        w_sub_all = [];

        for iCond = 1:nCond
            c = conds(iCond);

            if c > size(CMall,1) || b > size(CMall,2) || s > size(CMall,3)
                warning('Subject %s does not contain cond/band/side [%d/%d/%d]. Skipping.', subName, c, b, s);
                subValid = false; break
            end

            CM = CMall(c, b, s);

            if ~isfield(CM, 'rval') || isempty(CM.rval)
                warning('Missing rval for %s cond %d. Skipping.', subName, c);
                subValid = false; break
            end
            if ~isfield(CM, 'w') || isempty(CM.w)
                warning('Missing w for %s cond %d. Skipping.', subName, c);
                subValid = false; break
            end
            if ref > numel(CM.w) || isempty(CM.w{ref})
                warning('Missing w{%d} for %s cond %d. Skipping.', ref, subName, c);
                subValid = false; break
            end
            if ~isfield(CM, 't') || isempty(CM.t)
                warning('Missing t for %s cond %d. Skipping.', subName, c);
                subValid = false; break
            end

            t_sub = CM.t(:);
            w_sub = squeeze(CM.w{ref});          % [nTime x nSensors]
            r_sub = extract_rval(CM.rval, ref);  % [nSensors x 1]

            % Time consistency
            if isempty(t_ref)
                t_ref = t_sub;
            elseif numel(t_sub) ~= numel(t_ref) || any(abs(t_sub - t_ref) > 1e-12)
                warning('Time vector mismatch for %s. Skipping.', subName);
                subValid = false; break
            end

            % Sensor consistency
            if ~isempty(all_r)
                if numel(r_sub) ~= size(all_r,1) || size(w_sub,2) ~= size(all_w,2)
                    warning('Sensor count mismatch for %s. Skipping.', subName);
                    subValid = false; break
                end
            end

            r_sub_all(:, iCond) = r_sub;         % [nSensors x nCond]
            w_sub_all(:,:, iCond) = w_sub;        % [nTime x nSensors x nCond]
        end

        if ~subValid
            skippedCount = skippedCount + 1;
            continue
        end

        % Store subject
        loadedCount = loadedCount + 1;
        subjectIDs{loadedCount, 1} = subName;

        all_r(:, loadedCount, :) = r_sub_all;    % [nSensors x nSubjects x nCond]
        all_w(:,:, loadedCount, :) = w_sub_all;  % [nTime x nSensors x nSubjects x nCond]

    catch ME
        warning('Failed to load/process %s: %s', matPath, ME.message);
        skippedCount = skippedCount + 1;
    end
end % subjects

if loadedCount == 0
    error('No valid subjects could be loaded from: %s', rootDir);
end

% =========================================================================
% 3. Extract labels from loaded cfg
% =========================================================================
condLabels = strings(1, nCond);
for iCond = 1:nCond
    c = conds(iCond);
    condLabels(iCond) = sprintf('Cond %d', c);
    if ~isempty(cfg_loaded) && isfield(cfg_loaded, 'conditions') && numel(cfg_loaded.conditions) >= c
        condLabels(iCond) = string(cfg_loaded.conditions(c).label);
    end
end

bandLabel = sprintf('Band %d', b);
sideLabel = sprintf('Side %d', s);
refLabel  = sprintf('Ref %d',  ref);

if ~isempty(cfg_loaded)
    if isfield(cfg_loaded,'trf') && isfield(cfg_loaded.trf,'fwd_bands') && ...
       isfield(cfg_loaded.trf.fwd_bands,'label') && numel(cfg_loaded.trf.fwd_bands) >= b
        bandLabel = char(string(cfg_loaded.trf.fwd_bands(b).label));
    end
    if isfield(cfg_loaded,'trf') && isfield(cfg_loaded.trf,'sides') && numel(cfg_loaded.trf.sides) >= s
        sideLabel = char(string(cfg_loaded.trf.sides{s}));
    end
    if isfield(cfg_loaded,'ref_sources') && size(cfg_loaded.ref_sources,1) >= ref
        refLabel = char(string(cfg_loaded.ref_sources{ref,1}));
    end
end

% =========================================================================
% 4. Group averages  (one entry per condition)
% =========================================================================
% Preallocate
nTime    = size(all_w, 1);
nSensors = size(all_w, 2);
nSens_r  = size(CM.rval, 1);

r_group_mean = nan(nSens_r, nCond);   % [nSensors x nCond]
w_group_grand = nan(nTime, nCond);     % [nTime x nCond]
w_roi_mean    = nan(nTime, nCond);
w_roi_sem     = nan(nTime, nCond);
w_grand_sem   = nan(nTime, nCond);

% Sensor ROI
if isempty(cfg_plot.sensors)
    sensorIdx   = 1:nSensors;
    sensorLabel = 'All sensors';
else
    sensorIdx = cfg_plot.sensors(:)';
    sensorIdx = sensorIdx(sensorIdx >= 1 & sensorIdx <= nSensors);
    if isempty(sensorIdx)
        error('cfg_plot.sensors does not contain valid sensor indices.');
    end
    sensorLabel = sprintf('ROI (%d sensors)', numel(sensorIdx));
end

% Common zlim across conditions (computed after averaging)
r_all_conds = nan(nSens_r, nCond);

for iCond = 1:nCond

    r_c = squeeze(all_r(:,:,iCond));                      % [nSensors x nSubjects]
    w_c = squeeze(all_w(:,:,:,iCond));                    % [nTime x nSensors x nSubjects]

    r_group_mean(:, iCond) = mean(r_c, 2, 'omitnan');
    r_all_conds(:,  iCond) = r_group_mean(:, iCond);

    w_mean_c  = mean(w_c, 3, 'omitnan');                  % [nTime x nSensors]
    w_grand_c = mean(w_mean_c, 2, 'omitnan');             % [nTime x 1]

    w_sub_sensor = squeeze(mean(w_c, 2, 'omitnan'));      % [nTime x nSubjects]
    w_group_grand(:, iCond) = w_grand_c;
    w_grand_sem(:,   iCond) = std(w_sub_sensor, 0, 2, 'omitnan') ./ sqrt(loadedCount);

    w_roi_sub = squeeze(mean(w_c(:, sensorIdx, :), 2, 'omitnan'));  % [nTime x nSubjects]
    w_roi_mean(:, iCond) = mean(w_roi_sub, 2, 'omitnan');
    w_roi_sem(:,  iCond) = std(w_roi_sub,  0, 2, 'omitnan') ./ sqrt(loadedCount);
end

% Common zlim for topoplot
validAll = r_all_conds(~isnan(r_all_conds));
if isempty(validAll)
    zlim_common = [-1 1] * 1e-6;
else
    zlim_common = [min(validAll) max(validAll)];
    if diff(zlim_common) == 0
        zlim_common = zlim_common + [-1 1] * 1e-6;
    end
end

% =========================================================================
% 5. Figure  — subplot(nCond, 3, ...)
% =========================================================================
fig = figure( ...
    'Name',     'Group TRF', ...
    'Color',    'w', ...
    'Position', [100 100 1200 350*nCond]);

for iCond = 1:nCond

    rowOffset = (iCond - 1) * 3;

    % ---------------------------------------------------------------------
    % Panel 1 : Group mean TRF (all sensors)
    % ---------------------------------------------------------------------
    colors = lines(nCond);

    subplot(nCond, 3, rowOffset + 1);

    plot_with_sem(t_ref, w_group_grand(:,iCond), w_grand_sem(:,iCond), colors(iCond,:));
    xline(0, '--k', 'LineWidth', 1);
    xlabel('Time lag (ms)');
    ylabel('TRF weight');
    title(sprintf('Group TRF — %s', condLabels(iCond)), 'Interpreter', 'none');
    grid on; axis tight; box off;

    % ---------------------------------------------------------------------
    % Panel 2 : ROI TRF
    % ---------------------------------------------------------------------
    subplot(nCond, 3, rowOffset + 2);

    plot_with_sem(t_ref, w_roi_mean(:,iCond), w_roi_sem(:,iCond), colors(iCond,:));
    xline(0, '--k', 'LineWidth', 1);
    xlabel('Time lag (ms)');
    ylabel('TRF weight');
    title(sprintf('%s — %s', sensorLabel, condLabels(iCond)), 'Interpreter', 'none');
    grid on; axis tight; box off;

    % ---------------------------------------------------------------------
    % Panel 3 : Topoplot
    % ---------------------------------------------------------------------
    subplot(nCond, 3, rowOffset + 3);

    try
        cfg_plot.topoplot_fun(r_group_mean(:, iCond), zlim_common);
    catch ME
        warning('Topoplot failed for cond %d: %s\nUsing fallback bar plot.', iCond, ME.message);
        bar(r_group_mean(:, iCond));
        xlabel('Sensors');
        ylabel('r');
    end
    title(sprintf('Prediction r — %s', condLabels(iCond)), 'Interpreter', 'none');
    colorbar;
    box off;

end % conditions

% -------------------------------------------------------------------------
% Global title
% -------------------------------------------------------------------------
sgtitle(sprintf('Group TRF | N=%d | Band: %s | Side: %s | Ref: %s', ...
    loadedCount, bandLabel, sideLabel, refLabel), ...
    'FontWeight', 'bold', 'Interpreter', 'none');

% -------------------------------------------------------------------------
% Annotation box
% -------------------------------------------------------------------------
annotationText = sprintf(['Subjects loaded : %d\n' ...
                          'Subjects skipped: %d\n' ...
                          'Mean r (all cond): %.4f\n' ...
                          'Max  r (all cond): %.4f'], ...
                          loadedCount, skippedCount, ...
                          mean(r_all_conds(:), 'omitnan'), ...
                          max(r_all_conds(:)));

annotation(fig, 'textbox', [0.86 0.88 0.13 0.10], ...
    'String',          annotationText, ...
    'FitBoxToText',    'on', ...
    'BackgroundColor', 'w', ...
    'EdgeColor',       [0.8 0.8 0.8], ...
    'FontSize',        9);

% =========================================================================
% 6. Console output
% =========================================================================
if cfg_plot.show_console
    fprintf('\n============================================================\n');
    fprintf('GROUP TRF SUMMARY\n');
    fprintf('------------------------------------------------------------\n');
    fprintf('Root directory : %s\n', rootDir);
    fprintf('Band           : %s\n', bandLabel);
    fprintf('Side           : %s\n', sideLabel);
    fprintf('Reference      : %s\n', refLabel);
    fprintf('Subjects used  : %d\n', loadedCount);
    fprintf('Subjects skip  : %d\n', skippedCount);
    for iCond = 1:nCond
        fprintf('--- %s ---\n', condLabels(iCond));
        fprintf('  Mean r : %.4f\n', mean(r_group_mean(:,iCond), 'omitnan'));
        fprintf('  Max  r : %.4f\n', max(r_group_mean(:,iCond)));
    end
    fprintf('============================================================\n');
end

% =========================================================================
% 7. Save figure
% =========================================================================
if cfg_plot.save
    if ~exist(cfg_plot.save_dir, 'dir')
        mkdir(cfg_plot.save_dir);
    end
    fileBase = sprintf('Group_TRF_B%d_S%d_R%d', b, s, ref);
    saveas(fig, fullfile(cfg_plot.save_dir, [fileBase '.png']));
    saveas(fig, fullfile(cfg_plot.save_dir, [fileBase '.fig']));
end

end % plot_group_TRF()

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function r = extract_rval(rval, ref)
sz = size(rval);
if ndims(rval) ~= 3
    error('Expected rval to be 3D [nSensors x nRefs x nFolds], got size %s.', mat2str(sz));
end
if ref > sz(2)
    error('Requested ref (%d) exceeds available references (%d).', ref, sz(2));
end
r = squeeze(rval(:, ref, :));
r = mean(r, 2, 'omitnan');
end



function plot_with_sem(x, y, sem, color)
if nargin < 4 || isempty(color)
    color = [1 0 0];   % rouge par défaut (votre choix actuel)
end
x = x(:); y = y(:); sem = sem(:);

fill([x; flipud(x)], [y-sem; flipud(y+sem)], ...
     min(color * 1.5 + 0.4, 1), ...   % version éclaircie pour la SEM
     'EdgeColor', 'none', 'FaceAlpha', 0.6);
hold on;
plot(x, y, 'Color', color, 'LineWidth', 2);
end