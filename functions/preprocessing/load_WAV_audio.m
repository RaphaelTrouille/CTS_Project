function WAV_audio = load_WAV_audio(audio_file, target_fs)

    WAV_audio = struct('signal', [], 'Fs', []);

    if exist('wavread', 'file')
        [sig, fs_orig, ~] = wavread(audio_file);
    else
        [sig, fs_orig] = audioread(audio_file);
    end

    sig = sig(:)';

    if nargin > 1 && ~isempty(target_fs) && fs_orig ~= target_fs
        sig =resample(sig, target_fs, fs_orig);
        WAV_audio.Fs = target_fs;
    else
        WAV_audio.Fs = fs_orig;
    end
    WAV_audio.signal = sig;
end