function center_frequencies = get_mel_freqs(lower_freq, upper_freq, nb_bands)
% GET_MEL_FREQS  Generate MEL-spaced center frequencies between two bounds.
%
% DESCRIPTION:
%   Converts the lower and upper frequency bounds to the MEL scale,
%   creates nb_bands linearly spaced values in MEL space, then converts
%   them back to Hz. This produces perceptually uniform frequency bands
%   that mimic the non-linear frequency resolution of the human auditory system.
%
% INPUTS:
%   lower_freq  - Lower frequency bound in Hz (e.g. 100)
%   upper_freq  - Upper frequency bound in Hz (e.g. 8000)
%   nb_bands    - Number of frequency bands / center frequencies to generate
%
% OUTPUTS:
%   center_frequencies  - Vector of nb_bands center frequencies in Hz,
%                         linearly spaced on the MEL scale
%
% FORMULA:
%   MEL = 2595 * log10(1 + f / 700)
%   f   = 700 * (10^(MEL / 2595) - 1)
%
% USAGE:
%   cf = get_mel_freqs(100, 8000, 28);
%
% -------------------------------------------------------------------------    
    
    % Convert bounds to MEL scale and interpolate linearity
    center_frequencies = [lower_freq upper_freq];
    Fmel = 2595 * log10(1 + center_frequencies / 700);
    Fmel = linspace(Fmel(1), Fmel(2), nb_bands);
    
    % Convert back to Hz
    center_frequencies = 700 * (10.^(Fmel/2595)-1);
end