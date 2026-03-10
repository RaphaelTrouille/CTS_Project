function [t1, t2, L, gof] = alignement_verification(tds, dec, MISC_struct, audio_struct, perm_stat, do_plot)
% ALIGNEMENT_VERIFICATION    Verify the temporal alignment between a MEG signal
%                            and a downsampled audio signal.
%
% DESCRIPTION:
%   This function checks the co-registration quality between a MEG-rate signal
%   (MISC_struct) and an audio signal (audio_struct) by comparing their
%   bandpass-filtered versions in the 50-330 Hz bandwidth
%   on each respective timeline. A temporal offset (dec) is applied to account
%   for a known or estimated delay. Alignment quality is assessed via Pearson
%   correlation on the absolute values of both signals. Optionally plots both
%   signals for visual inspection.
%
% INPUTS:
%   tds          - Downsampling indices into the audio signal (mapping audio
%                  samples to MEG timepoints)
%   dec          - Temporal offset in samples (positive = MEG leads audio)
%   MISC_struct  - Structure with fields:
%                    .signal : MEG-rate signal vector
%                    .Fs     : MEG sampling rate in Hz
%   audio_struct - Structure with fields:
%                    .signal : audio signal vector
%                    .Fs     : audio sampling rate in Hz
%   perm_stat    - Boolean flag; if true, skips correlation computation
%                  (used during permutation testing)
%   do_plot      - Boolean flag; if true, plots both signals for visual check
%
% OUTPUTS:
%   t1   - Indices into MISC_struct.signal for the aligned segment
%   t2   - Indices into the downsampled audio signal for the aligned segment
%   L    - Length of the aligned overlap segment (in samples)
%   gof  - Goodness-of-fit: Pearson correlation between |MISCf(t1)| and
%          |Yds(t2)|. Returns [] if perm_stat is true. Should be above 0.5.
%
% METHOD:
%   Both signals are bandpass filtered between 50-330 Hz (cosine filter)
%   on their respective sampling rates before comparison. The audio signal
%   is then downsampled to MEG rate via tds indices. The overlap segment
%   of length L is determined by dec and signal lengths.
%
% DEPENDENCIES:
%   - cosine_filter.m (custom function, must be on MATLAB path)
%
% USAGE:
%   [t1, t2, L, gof] = alignement_verification(tds, 0, MISC_struct, audio_struct, false, true);
%
% -------------------------------------------------------------------------    
    % Discard downsampling indices that exceed the audio signal length
    tds = tds(tds<length(audio_struct.signal));

    % Bandpass filter the MEG-rate signal (50 - 330 Hz)
    F_h_x = cosine_filter(length(MISC_struct.signal),...
                         {'high' 'low'},...
                         [50 330]/MISC_struct.Fs*2,...
                         [5 5]/MISC_struct.Fs*2);
    MISCf = real(...
                ifft(...
                    fft(...
                        MISC_struct.signal).*F_h_x)...
                        );

    % Bandpass filter the audio signal (50 - 330 Hz)
    F_h_x = cosine_filter(length(audio_struct.signal),...
                          {'high' 'low'}, ...
                          [50 330] / audio_struct.Fs*2, ...
                          [5 5] / audio_struct.Fs*2);
    
    Yf = real(...
            ifft(...
                fft( ...
                    audio_struct.signal).*F_h_x)...
                    );

    % Downsample filtered audio to MEG rate
    Yds = Yf(tds);
    
    % Compute aligned overlap segment accounting for temporal offset dec
    L = min(length(MISC_struct.signal) - max(dec, 0),...
            length(Yds) + max(-dec, 0));
    t1 = max(dec, 0)+(1:L);
    t2 = max(-dec, 0)+(1:L);

    % Optional: plot both signals for visual alignment check
    if do_plot == true
        plot(MISCf(t1));    hold on
        plot(Yds(t2)/norm(Yds(t2)) * norm(MISCf(t1)), 'g');
    end

    % Compute goodness-of-fit (skipped during permutation testing)
    if ~perm_stat
        gof = corr(abs(MISCf(t1))', abs(Yds(t2))');
        disp(['g.o.f of the sound coregistration = ', num2str(gof), '. Should be above 0.5'])
    end

end