function [dec, tds] = realign_sound_file(syncfile, MISCaudio, WAVaudio)

    if isfile(syncfile)
        data = load(syncfile, 'dec', 'tds');
        dec = data.dec;
        tds = data.tds;
    else
        [dec, tds] = CM_CVC_realign_MISC_son(MISCaudio.signal,...
                                            WAVaudio.signal,...
                                            MISCaudio.Fs,...
                                            WAVaudio.Fs);
        save(syncfile, 'dec', 'tds');
    end
end