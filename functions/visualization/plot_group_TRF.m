function plot_group_TRF(rootDir, cfg_plot)

% =========================================================================
% 0) Input checks and defaults
% =========================================================================
validateattributes(rootDir, {'char', 'string'}, {'nonempty'}, mfilename, 'rootDir', 1);
root = char(rootDir);

validateattributes(cfg_plot, {'struct'}, {'nonempty'}, mfilename, 'cfg_plot', 2);

requiredFields = {'cond', 'band', 'ref'};
for iF = 1:numel(requiredFields)
    if ~isfield(cfg_plot, requiredFields{iF})
        error('Missing required cfg_plot field: cfg_plot.%s', requiredFields{iF});
    end
end

validateattributes(cfg_plot.cond, {'numeric'}, {'scalar', 'integer', 'positive'});
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
    cfg_plot.pattern = '*.mat';
end
if ~isfield(cfg_plot, 'sensors') || isempty(cfg_plot.sensors)
    cfg_plot.sensors = [];
end
if ~isfield(cfg_plot, 'topoplot_fun') || isempty(cfg_plot.topoplot_fun)
    cfg_plot.topoplot_fun = @FT_quick_topoplot;
end

c = cfg_plot.cond;
b = cfg_plot.band;
s = cfg_plot.side;
ref = cfg_plot.ref;

if ~exist(rootDir, 'dir')
    error('Input directory does not exist: %s', rootDir);
end

% =========================================================================
% 1. Find subject folder and result files
% =========================================================================

subDirs = dir(rootDir);
subDirs = subDirs([subDirs.isdir]);
subDirs = subDirs(~ismember({subDirs.name}, {'.', '..'}));

if isempty(subDirs)
    error('No subject folders found in: %s', rootDir);
end

subjectIDs   = {};
all_r        = [];
all_w        = [];
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
        warning('No .mat file found in %s (pattern: %s). Skippiing.', subPath, cfg_plot.pattern)
        skippedCount = skippedCount + 1;
        continue
    end

    % Use first matching file by default
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

        % Safety on dimensions
        if c > size(CMall, 1) || b > size(CMall, 2) || s > size(CMall, 3)
            warning('Subjects %s does not contain requested cond/band/side. Skipping.', subName);
            skippedCount = skippedCount + 1;
            continue
        end

        CM = CMall(c,b,s);

        % -----------------------------------------------------------------
        % Extract r
        % -----------------------------------------------------------------
        if ~isfield(CM, 'rval') || isempty(CM.rval)
            warning('Missing rval fo %s. Skipping.', subName);
            skippedCount = skippedCount + 1;
            continue
        end

        r_sub = extract_rval(CM.rval, ref);      % [nSensors x 1]
        
        % -----------------------------------------------------------------
        % Extract w
        % -----------------------------------------------------------------
        if ~isfield(CM, 'w') || isempty(CM.w)
            warning('Missing w fo %s. Skipping.', subName);
            skippedCount = skippedCount + 1;
            continue
        end

        if ref > numel(CM.w) || isempty(CM.w{ref})
            warning('Missing w{%d} for %s. Skipping.', ref, subName);
            skippedCount = skippedCount + 1;
            continue
        end

        if ~isfield(CM, 't') || isempty(CM.t)
            warning('Missing t fo %s. Skipping.', subName);
            skippedCount = skippedCount + 1;
            continue
        end

        t_sub = CM.t(:);
        w_sub = squeeze(CM.w{ref});     % [nTime x nSensors]

        % -----------------------------------------------------------------
        % Time consistency
        % -----------------------------------------------------------------
        if isempty(t_ref)
            t_ref = t_sub;
        else
            if numel(t_sub) ~= numel(t_ref) || any(abs(t_sub - t_ref) > 1e-12)
                warning('Time vector mismatch for %s. Skipping.', subName);
                skippedCount = skippedCount + 1;
                continue
            end
        end

        % -----------------------------------------------------------------
        % Sensor consistency
        % -----------------------------------------------------------------
        if isempty(all_r)
            nSensors_ref = numel(r_sub);
        else
            if numel(r_sub) ~= size(all_r,1) || size(w_sub, 2) ~= size(all_w, 2)
                warning('Sensor count mismatch for %s. Skipping.', subName);
                skippedCount = skippedCount + 1;
                continue
            end
        end

        % Store
        loadedCount = loadedCount + 1;
        subjectIDs{loadedCount, 1} = subName;

        all_r(:, loadedCount)   = r_sub;
        all_w(:,:, loadedCount) = w_sub;
    
    catch ME
        warning('Failed to load/process %s: %s', matPath, ME.message);
        skippedCount = skippedCount + 1;
        continue
    end % Try
end % Subjects 


if loadedCount == 0
    error('No valid subjects could be loaded from: %s', rootDir)
end

% =========================================================================
% 3) Extract labels from loaded cfg structure
% =========================================================================

condLabel = sprintf('Cond %d', c);
bandLabel = sprintf('Band %d', b);
sideLabel = sprintf('Side %d', s);
refLabel  = sprintf('Ref %d', ref);

if ~isempty(cfg_loaded)

    % -------------------------
    % Condition labels
    % -------------------------
    if isfield(cfg_loaded, 'conditions') && ...
       isfield(cfg_loaded.conditions, 'label') && ...
       numel(cfg_loaded.conditions.label) >= c

        condLabel = char(string(cfg_loaded.conditions(c).label));
    end

    % -------------------------
    % Band labels
    % -------------------------
    if isfield(cfg_loaded, 'trf') && ...
       isfield(cfg_loaded.trf, 'fwd_bands') && ...
       isfield(cfg_loaded.trf.fwd_bands, 'label') && ...
       numel(cfg_loaded.trf.fwd_bands.label) >= b

        bandLabel = char(string(cfg_loaded.trf.fwd_bands(b).label));
    end

    % -------------------------
    % Side labels
    % -------------------------
    if isfield(cfg_loaded, 'trf') && ...
       isfield(cfg_loaded.trf, 'sides') && ...
       numel(cfg_loaded.trf.sides) >= s

        sideLabel = char(string(cfg_loaded.trf.sides{s}));
    end

    % -------------------------
    % Reference labels
    % -------------------------
    if isfield(cfg_loaded, 'ref_sources') && ...
       size(cfg_loaded.ref_sources,1) >= ref

        refLabel = char(string(cfg_loaded.ref_sources{ref,1}));
    end

end

% =========================================================================
% 4. Group averages
% =========================================================================

% Group topography
r_group_mean = mean(all_r, 2, 'omitnan');
r_group_sem  = std(all_r, 0, 2, 'omitnan') ./sqrt(loadedCount);

% Group TRF
w_group_mean  = mean(all_w, 3, 'omitnan');       % [nTime x nSensors]
w_group_grand = mean(w_group_mean, 2, 'omitnan'); % [nTime x 1] 

% Subject mean across sensors
w_sub_sensorMean = squeeze(mean(all_w, 2, 'omitnan')); % [nTime x nSubjects]

% ROI or all-sensor average
if isempty(cfg_plot.sensors)
    sensorIdx   = 1:size(w_group_mean, 2);
    sensorLabel = 'All sensors';
else
    sensorIdx = cfg_plot.sensors(:)';
    sensorIdx = sensorIdx(sensorIdx >= 1 & sensorIdx <= size(w_group_mean, 2));
    if isempty(sensorIdx)
        error('cfg_plot.sensors does not contain valid sensor indices.');
    end
    sensorLabel = sprintf('ROI (%d sensors)', numel(sensorIdx));

end

w_roi_subjects = squeeze(mean(all_w(:, sensorIdx, :), 2, 'omitnan')); % [nTime x nSubjects]
w_roi_mean = mean(w_roi_subjects, 2, 'omitnan');
w_roi_sem  = std(w_roi_subjects, 0, 2, 'omitnan')./ sqrt(loadedCount);

% =========================================================================
% 5. Figure
% =========================================================================

fig = figure( ...
    'Name', 'Group TRF', ...
    'Color', 'w', ...
    'Position', [100 100 1450 450]);

% -------------------------------------------------------------------------
% Panel 1: Group mean TRF across all sensors
% -------------------------------------------------------------------------
subplot(1, 3, 1)

plot_with_sem(t_ref, w_group_grand, std(w_sub_sensorMean, 0, 2, 'omitnan') ./ sqrt(loadedCount));
hold on;
xline(0, '--k', 'LineWidth', 1);

xlabel('Time lag (ms)');
ylabel('TRF weight');
title(sprintf('Group mean TRF\nCond: %s | Band: %s\nSide: %s | Ref: %s', ...
    condLabel, bandLabel, sideLabel, refLabel), ...
    'Interpreter', 'none');
grid on;
axis tight;
box off;

% -------------------------------------------------------------------------
% Panel 2: ROI / selected TRF
% -------------------------------------------------------------------------
subplot(1, 3, 2);

plot_with_sem(t_ref, w_roi_mean, w_roi_sem);
hold on;
xline(0, '--k', 'LineWidth', 1);

xlabel('Time lag (ms)');
ylabel('TRF weight');
title(sprintf('Sensor summary | %s', sensorLabel), 'Interpreter', 'none');
grid on;
axis tight;
box off;

% -------------------------------------------------------------------------
% Panel 3: Group topography
% -------------------------------------------------------------------------
subplot(1, 3, 3);

validVals = r_group_mean(~isnan(r_group_mean));
if isempty(validVals)
    zlim = [-1 1] * 1e-6;
else
    zlim = [min(validVals) max(validVals)];
    if diff(zlim ==0)
        zlim = zlim + [-1 1] * 1e-6;
    end
end

try 
    cfg_plot.topoplot_fun(r_group_mean, zlim);
catch ME
warning('Topoplot failed at group level: %s\nUsing fallback bar plot.', ME.message);
bar(r_group_mean);
xlabel('Sensors')
ylabel('Prediction performance (r)');
end

title(sprintf('Group prediction performance | %s', refLabel), ...
    'Interpreter','none');
colorbar;
box off;

% -------------------------------------------------------------------------
% Global title
% -------------------------------------------------------------------------
sgtitle(sprintf('Group TRF summary | N = %d subjects', loadedCount), ...
        'FontWeight', 'bold', ...
        'Interpreter', 'none');

annotationText = sprintf(['Subjects loaded : %d\n' ...
                          'Subjects skipped: %d\n' ...
                          'Mean r          : %.4f\n' ...
                          'Max r           : %.4f'], ...
                          loadedCount, skippedCount, ...
                          mean(r_group_mean, 'omitnan'), ...
                          max(r_group_mean));

annotation(fig, 'textbox', [0.86 0.68 0.13 0.18], ...
           'String', annotationText, ...
           'FitBoxToText', 'on', ...
           'BackgroundColor', 'w', ...
           'EdgeColor', [0.8 0.8 0.8], ...
           'FontSize', 9);

% =========================================================================
% 6. Console output
% =========================================================================

if cfg_plot.show_console
    fprintf('\n============================================================\n');
    fprintf('GROUP TRF SUMMARY\n');
    fprintf('------------------------------------------------------------\n');
    fprintf('Root directory : %s\n', rootDir);
    fprintf('Condition      : %s\n', condLabel);
    fprintf('Band           : %s\n', bandLabel);
    fprintf('Side           : %s\n', sideLabel);
    fprintf('Reference      : %s\n', refLabel);
    fprintf('Subjects used  : %d\n', loadedCount);
    fprintf('Subjects skip  : %d\n', skippedCount);
    fprintf('Mean group r   : %.4f\n', mean(r_group_mean, 'omitnan'));
    fprintf('Max group r    : %.4f\n', max(r_group_mean));
    fprintf('============================================================\n');
end

% =========================================================================
% 7. Save figure
% =========================================================================

if cfg_plot.save
    if ~exist(cfg_plot.save_dir, 'dir')
        mkdir(cfg_plot.save_dir);
    end

    fileBase = sprintf('Group_TRF_C%d_B%d_S%d_R%d', c, b, s, ref);

    saveas(fig, fullfile(cfg_plot.save_dir, [fileBase '.png']));
    saveas(fig, fullfile(cfg_plot.save_dir, [fileBase '.fig']));
end

end % plot_group_TRF()

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function r = extract_rval(rval, ref)
%EXTRACT_RVAL Extract mean prediction performance across folds.
%
%   INPUT:
%       rval : [nSensors x nRefs x nFolds]
%
%   OUTPUT:
%       r    : [nSensors x 1]

sz = size(rval);

if ndims(rval) ~= 3
    error('Expected rval to be 3D [nSensors x nRefs x nFolds], got size %s.', mat2str(sz));
end

nRefs = sz(2);
if ref > nRefs
    error('Requested ref (%d) exceeds available references (%d).', ref, nRefs);
end

r = squeeze(rval(:, ref, :));   % [nSensors x nFolds]
r = mean(r, 2, 'omitnan');      % [nSensors x 1]

end % extract_rval()


function plot_with_sem(x, y, sem)
%PLOT_WITH_SEM Plot mean curve with shaded SEM.

x = x(:);
y = y(:);
sem = sem(:);

fill([x; flipud(x)], [y-sem; flipud(y+sem)], ...
     [0.85 0.85 0.85], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.6);
hold on;
plot(x, y, 'k', 'LineWidth', 2);

end

