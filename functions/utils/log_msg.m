function log_msg(cfg, varargin)
% LOG_MSG  Conditional logging utility for the CTS pipeline.
%
% DESCRIPTION:
%   Prints a formatted message to the console only if cfg.verbose is true.
%   Supports the same syntax as fprintf (format string + arguments).
%   Optionally writes to a log file if cfg.log.file is defined.
%
% INPUTS:
%   cfg       - Pipeline configuration struct with fields:
%                 .verbose       : boolean, print to console if true
%                 .log.do        : (optional) boolean, write to file if true
%                 .log.file      : (optional) full path to the log file
%   varargin  - Format string and arguments, passed directly to fprintf
%               (same syntax as fprintf, e.g. '%s done in %.1f sec\n', name, t)
%
% USAGE:
%   log_msg(cfg, '[%d/%d] Processing: %s\n', n_sub, n_total, sub_name);
%   log_msg(cfg, '[WARNING] File not found: %s\n', filepath);
%
% -------------------------------------------------------------------------

    % Print to console
    if cfg.verbose
        fprintf(varargin{:});
    end

    % Write to log file if configured
    if isfield(cfg, 'log') && isfield(cfg.log, 'do') && cfg.log.do
        fid = fopen(cfg.log.file, 'a');  % Append mode
        if fid == -1
            warning('log_msg: Could not open log file: %s', cfg.log.file);
            return
        end
        % Prepend timestamp to each log entry
        fprintf(fid, '[%s] ', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        fprintf(fid, varargin{:});
        fclose(fid);
    end
end