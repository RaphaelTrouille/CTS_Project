function cfg_plot = plot_config()
cfg_plot.cond = 1;
cfg_plot.band = 1;
cfg_plot.side = 3;
cfg_plot.ref  = 2;

cfg_plot.save = true;
cfg_plot.save_dir = fullfile(pwd, 'results', 'figures');

cfg_plot.band_labels = {'delta','theta'};
cfg_plot.side_labels = {'left','right','both'};
cfg_plot.ref_labels  = {'Global','Attended','Noisy', 'Mouth aperture'};

cfg_plot.topoplot_fun = @FT_quick_topoplot;
cfg_plot.show_console = true;

end