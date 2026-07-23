function permuted_signal = permute_signal_segments(signal, Fs)
% PERMUTE_SIGNAL_SEGMENTS  Permute a speech signal by swapping its two
%                           active halves around a silent midpoint.
%
% DESCRIPTION:
%   This function performs a structured permutation of a speech signal
%   for use in permutation testing. It identifies the active speech region
%   (above a low-amplitude threshold), finds a natural silence point near
%   the temporal midpoint of that region, and swaps the two halves of the
%   speech segment. The silent margins before and after the speech are
%   preserved in place. This creates a temporally shuffled version of the
%   signal that preserves its spectral properties while disrupting temporal
%   alignment with neural responses.
%
% INPUTS:
%   signal  - Input signal vector (e.g. auditory envelope), row or column
%   Fs      - Sampling rate in Hz (used to define the search window around
%             the midpoint: +/- 10 seconds)
%
% OUTPUTS:
%   permuted_signal  - Signal with the two speech halves swapped, same
%                      length as input
%
% METHOD:
%   1. Detect speech onset (talk) and offset (tlop) as first/last samples
%      exceeding threshold (mean/20)
%   2. Find the temporal midpoint between onset and offset
%   3. Within a +/- 10s window around the midpoint, locate the local
%      minimum (natural silence) to use as the split point
%   4. Swap the two halves: [silence | half2 | half1 | silence]
%
% USAGE:
%   perm_env = permute_signal_segments(envelope, 1000);
%
% -------------------------------------------------------------------------
 
    sig = signal;

    % Detect active speech region using a low amplitude threshold (mean/20)
    th = mean(sig)/20;
    talk = find(sig>th, 1);         % Speech onset sample
    tlop = find(sig>th, 1, 'last'); % Speech offset sample
    
    % Find temporal midpoint between onset and offset
    tmid = round((talk + tlop)/2);

    % Refine midpoint to the local minimum within a +/-10s window (natural silence)
    [~, tcorr] = min(sig(tmid-10*Fs:tmid+10*Fs));
    tmid = tmid + tcorr - 10 * Fs-1;

    % Build permuted index vector: swap the two speech halves
    t_shuffle = (1:length(sig));
    
    t_shuffle = [t_shuffle(1:talk-1) ...      % Pre-speech silence (unchanged)
                 t_shuffle(tmid:tlop - 1) ... % Second half of speech (moved first) 
                 t_shuffle(talk:tmid-1) ...   % Fisrt half of speech (moved second)
                 t_shuffle(tlop:end)];        % Post-speech silence (unchanged)
    
    permuted_signal = signal(t_shuffle);
end