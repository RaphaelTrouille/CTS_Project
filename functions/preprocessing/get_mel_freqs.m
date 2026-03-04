function center_frequencies = get_mel_freqs(lower_freq, upper_freq, nb_bands)
    center_frequencies = [lower_freq upper_freq];
    Fmel = 2595 * log10(1 + center_frequencies / 700);
    Fmel = linspace(Fmel(1), Fmel(2), nb_bands);
    center_frequencies = 700 * (10.^(Fmel/2595)-1);
end