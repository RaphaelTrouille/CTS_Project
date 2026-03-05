function enveloppe = enveloppe_extraction(audio_struct, center_frequencies)
% ENVELOPPE_EXTRACTION  Compute the auditory envelope of a speech signal
%                       using a gammatone filterbank.
%
% DESCRIPTION:
%   This function extracts the auditory envelope of a speech signal by
%   filtering it through a bank of gammatone filters centered at specified
%   frequencies. For each frequency channel, the filter response is
%   amplitude-normalized, the half-wave rectified envelope is smoothed
%   with a Hanning kernel, and compressed with a power-law (exponent 0.6)
%   to approximate auditory nerve firing rates. Envelopes across all
%   channels are summed to produce the final broadband auditory envelope.
%
% INPUTS:
%   audio_struct        - Structure with fields:
%                           .signal : audio signal vector (1 x N)
%                           .Fs     : sampling rate in Hz
%   center_frequencies  - Vector of center frequencies (in Hz) for the
%                         gammatone filterbank (e.g. ERB-spaced frequencies)
%
% OUTPUTS:
%   enveloppe  - Broadband auditory envelope (1 x N), same length as
%                audio_struct.signal
%
% METHOD:
%   For each center frequency:
%     1. Design a gammatone filter (via AMT_toolbox_gammatone)
%     2. Normalize filter amplitude using its impulse response FFT magnitude
%     3. Filter the signal and take the absolute value (full-wave rectification)
%     4. Smooth with a Hanning kernel (window ~ Fs/200 samples, i.e. ~5 ms)
%     5. Apply power-law compression (x^0.6)
%     6. Accumulate across frequency channels
%
% DEPENDENCIES:
%   - AMT_toolbox_gammatone.m (AMToolbox wrapper, must be on MATLAB path)
%
% USAGE:
%   cf = get_mel_freqs(150, 7000, 31);   % ERB-spaced center frequencies
%   env = enveloppe_extraction(audio_struct, cf);
%

% -------------------------------------------------------------------------    
    
    Fs = audio_struct.Fs;
    signal = audio_struct.signal;
    enveloppe = zeros(size(signal));
    
    for freq = 1:length(center_frequencies)
        % Design gannatone filter for this frequency channel
        [b, a] = AMT_toolbox_gammatone(center_frequencies(freq), Fs);
        
        % Estimate filter gain via its impulse response in the frequency
        % domain
        test = zeros(1, Fs);
        test(1) = 1;
        test = filter(b, a, test);
        ampl = norm(abs(fft(test)))/100;
        
        % Filter the signal and normalize by estimated gain
        Yf = filter(b, a, signal)/ampl;
        
        % Smooth the rectified enveloppe with a Hanning window
        kern = hanning(round(Fs/200));
        kern = kern/sum(kern);
        Yf = conv(abs(Yf), kern, 'same');

        % Power-law compression (0.6 exponent) + accumulate across channels
        enveloppe = enveloppe + Yf.^0.6;
    end
end