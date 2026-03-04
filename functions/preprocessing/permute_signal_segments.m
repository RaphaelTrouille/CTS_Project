function permuted_signal = permute_signal_segments(signal, Fs)
    sig = signal;
    th = mean(sig)/20;
    talk = find(sig>th, 1);
    tlop = find(sig>th, 1, 'last');
    tmid = round((talk + tlop)/2);
    [~, tcorr] = min(sig(tmid-10*Fs:tmid+10*Fs));
    tmid = tmid + tcorr - 10 * Fs-1;
    t_shuffle = (1:length(sig));
    t_shuffle = [t_shuffle(1:talk-1) t_shuffle(tmid:tlop - 1) t_shuffle(talk:tmid-1) t_shuffle(tlop:end)];
    permuted_signal = signal(t_shuffle);
end