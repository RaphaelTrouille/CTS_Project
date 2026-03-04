function enveloppe = enveloppe_extraction(audio_struct, center_frequencies)
    
    Fs = audio_struct.Fs;
    signal = audio_struct.signal;
    enveloppe = zeros(size(signal));
    
    for freq = 1:length(center_frequencies)
        [b, a] = AMT_toolbox_gammatone(center_frequencies(freq), Fs);
        test = zeros(1, Fs);
        test(1) = 1;
        test = filter(b, a, test);
        ampl = norm(abs(fft(test)))/100;
        Yf = filter(b, a, signal)/ampl;
        kern = hanning(round(Fs/200));
        kern = kern/sum(kern);
        Yf = conv(abs(Yf), kern, 'same');
        enveloppe = enveloppe + Yf.^0.6;
    end
end