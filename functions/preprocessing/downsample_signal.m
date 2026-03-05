function MISC = downsample_signal(envelope, tds, t1, t2, target_size)
% DOWNSAMPLE_SIGNAL  Downsample a signal and insert it into a target array.
%
% DESCRIPTION:
%   Extracts samples from an envelope signal at specified indices (tds),
%   then places them into a zero-initialized output array of a given size
%   at target indices (t1), using source indices (t2) for alignment.
%   This allows flexible downsampling with index-based mapping between
%   source and target timelines.
%
% INPUTS:
%   enveloppe    - Input signal vector to downsample
%   tds          - Indices into enveloppe to extract (downsampling positions)
%   t1           - Indices in the output array where samples will be inserted
%   t2           - Indices into the downsampled signal Yds to use as source
%   target_size  - Size reference for the output array (e.g. size of MEG signal)
%
% OUTPUTS:
%   MISC  - Zero-initialized array of size target_size, with downsampled
%           values inserted at positions t1
%
% USAGE:
%   MISC = downsample_signal(enveloppe, tds, t1, t2, size(meg_signal));
%
% -------------------------------------------------------------------------
    
    % Extract downsampled values from the envelope
    Yds = envelope(tds);
    
    % Insert into a zero-padded output array using index mapping
    MISC = zeros(size(target_size));
    MISC(t1) = Yds(t2);
end