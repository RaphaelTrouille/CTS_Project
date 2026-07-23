function WAV_audio = load_WAV_audio(audio_file, target_fs)
% LOAD_WAV_AUDIO  Load a WAV audio file and optionally resample it.
%
% DESCRIPTION:
%   Reads a WAV audio file and returns it as a standardized structure
%   containing the signal and its sampling rate. If a target sampling rate
%   is provided and differs from the original, the signal is resampled
%   accordingly. Handles both legacy MATLAB (wavread) and modern (audioread)
%   APIs for compatibility.
%
% INPUTS:
%   audio_file  - Full path to the WAV file (string)
%   target_fs   - (Optional) Target sampling rate in Hz. If omitted, empty,
%                 or equal to the file's native rate, no resampling is applied.
%
% OUTPUTS:
%   WAV_audio   - Structure with fields:
%                   .signal : audio signal as a row vector (1 x N)
%                   .Fs     : sampling rate in Hz (original or resampled)
%
% NOTES:
%   - The signal is always returned as a row vector regardless of the
%     original file format.
%   - wavread is used if available (MATLAB < R2015a), otherwise audioread
%     is used (MATLAB >= R2012b).
%
% USAGE:
%   audio = load_WAV_audio('path/to/file.wav', 1000);
%
% -------------------------------------------------------------------------
    
    WAV_audio = struct('signal', [], 'Fs', []);
    
    % Use legacy wavread if available, otherwise use audioread (>= R2012b)
    if exist('wavread', 'file')
        [sig, fs_orig, ~] = wavread(audio_file);
    else
        [sig, fs_orig] = audioread(audio_file);
    end

    % Ensure signal is a row vector
    sig = sig(:)';
    
    % Resample if a target sampling rate is specifid and differs from
    % original
    if nargin > 1 && ~isempty(target_fs) && fs_orig ~= target_fs
        sig = resample(sig, target_fs, fs_orig);
        WAV_audio.Fs = target_fs;
    else
        WAV_audio.Fs = fs_orig;
    end
    WAV_audio.signal = sig;
end