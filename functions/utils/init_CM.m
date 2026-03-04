function CM = init_CM(Fs, quantum, window_type)

CM = struct();
CM.win = window_type;
CM.quantum = 2 * Fs * quantum; % seconds
CM.ref = [];
CM.Find = 1:CM.quantum/Fs*20+1;
end
