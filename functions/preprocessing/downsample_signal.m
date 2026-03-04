function MISC = downsample_signal(enveloppe, tds, t1, t2, target_size)
Yds = enveloppe(tds);
MISC = zeros(size(target_size));
MISC(t1) = Yds(t2);
end