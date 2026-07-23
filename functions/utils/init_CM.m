function CM = init_CM(Fs, quantum, window_type)
% INIT_CM  Initialize a CM (Cortical Mapping) structure.
%
% DESCRIPTION:
%   Creates and returns an empty CM structure with core parameters required
%   for cortical mapping analysis. The quantum defines the analysis window
%   length, and Find precomputes the frequency index vector used in spectral
%   computations (up to 20 Hz).
%
% INPUTS:
%   Fs           - Sampling rate in Hz
%   quantum      - Analysis window duration in seconds (half-window;
%                  stored as 2 * Fs * quantum samples)
%   window_type  - Window type for spectral analysis (e.g. 'hanning', 'hamming')
%
% OUTPUTS:
%   CM  - Initialized CM structure with fields:
%           .win    : window type string
%           .quantum: analysis window length in samples (2 * Fs * quantum)
%           .ref    : empty reference signal array (populated by add_CM_ref)
%           .Find   : frequency index vector for 0-20 Hz spectral bins
%
% USAGE:
%   CM = init_CM(1000, 0.5, 'hanning');
%
% NOTE:
%   CM.Find spans from 1 to (quantum/Fs*20 + 1), corresponding to frequency
%   bins from 0 Hz up to 20 Hz given the window length. This assumes that
%   frequency resolution is 1/quantum Hz per bin.
%
% -------------------------------------------------------------------------
    CM = struct();
    CM.win = window_type;
    CM.quantum = 2 * Fs * quantum;  % Analysis window length in samples
    CM.ref = [];                    % Reference signals (see add_CM_ref)
    CM.Find = 1:CM.quantum/Fs*20+1; % Frequency indices for 0-20Hz bins
end
