function [dec, tds] = realign_sound_file(syncfile, MISCaudio, WAVaudio)
% REALIGN_SOUND_FILE  Compute or load the temporal alignment between a MISC
%                     channel and a WAV audio file.
%
% DESCRIPTION:
%   Estimates the sample offset (dec) and downsampling indices (tds) needed
%   to align a WAV audio signal with a MEG MISC channel signal. If a
%   previously computed synchronization file exists, alignment parameters
%   are loaded directly from it to avoid redundant computation. Otherwise,
%   alignment is computed via CM_CVC_realign_MISC_son and the result is
%   cached to disk for future use.
%
% INPUTS:
%   syncfile   - Full path to the .mat cache file for storing/loading
%                alignment parameters
%   MISCaudio  - Structure with fields:
%                  .signal : MISC channel signal (row vector)
%                  .Fs     : MISC sampling rate in Hz
%   WAVaudio   - Structure with fields:
%                  .signal : WAV audio signal (row vector)
%                  .Fs     : WAV sampling rate in Hz
%
% OUTPUTS:
%   dec  - Sample offset between the MISC and WAV signals (integer).
%          Positive values indicate the MEG signal leads the audio.
%   tds  - Downsampling indices mapping WAV samples to MISC timepoints
%
% DEPENDENCIES:
%   - CM_CVC_realign_MISC_son.m (cartographie_motrice toolbox)
%
% USAGE:
%   [dec, tds] = realign_sound_file('path/to/sync.mat', MISCaudio, WAVaudio);
%
% -------------------------------------------------------------------------
    
    % Load cached alignment if exists
    if isfile(syncfile)
        data = load(syncfile, 'dec', 'tds');
        dec = data.dec;
        tds = data.tds;
    else
        % Compute alignment and cache result to disk
        [dec, tds] = CM_CVC_realign_MISC_son(MISCaudio.signal,...
                                            WAVaudio.signal,...
                                            MISCaudio.Fs,...
                                            WAVaudio.Fs);
        save(syncfile, 'dec', 'tds');
    end
end