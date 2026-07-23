clear; clc;

cfg = plot_config();

if isfield(cfg, 'subjects') && strcmp(cfg.subjects, 'all')    
    sub_dirs = dir(fullfile(cfg.deriv_dir, 'meg_*'));
else
    for i = 1:length(cfg.subjects)
        sub_dirs{i} = fullfile(cfg.deriv_dir, cfg.subjects{i});
    end
end

for i = 1:length(sub_dirs)
    sub_name = sub_dirs{i}.name;
    file_to_load = fullfile(cfg.deriv_dir, sub_name, [sub_name cfg.file_name]);
    if exist(file_to_load, 'file')
        fprintf('Processing visualization for %s...\n', sub_name)
        load(file_to_load, 'CMall');
        % Process the data as needed

        plot_TRF(CMall, cfg, sub_name)
    else
        warning('File %s does not exist.', file_to_load);
    end
end