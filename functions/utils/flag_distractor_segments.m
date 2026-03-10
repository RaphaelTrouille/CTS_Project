function bad = flag_distractor_segments(bad, t_dis, tds, dec, Fs)
% FLAG_DISTRACTOR_SEGMENTS  Mark distractor time segments as artefacts
%                            in the bad channel mask.
%
% DESCRIPTION:
%   For each distractor interval (e.g. a singing segment) defined in t_dis,
%   this function finds the closest corresponding samples in the MEG timeline
%   using the downsampling index vector (tds), then flags those samples as
%   bad in the artefact mask.
%
% INPUTS:
%   bad    - Binary artefact mask (1 x N); 1 = artefact, 0 = clean.
%            Updated in-place and returned.
%   t_dis  - M x 2 matrix of distractor intervals in seconds,
%            one row per distractor: [onset, offset]
%   tds    - Downsampling index vector mapping audio samples to MEG timepoints
%            (output of realign_sound_file)
%   dec    - Temporal offset in samples between MEG and audio timelines
%            (output of realign_sound_file)
%   Fs     - MEG sampling rate in Hz
%
% OUTPUTS:
%   bad    - Updated artefact mask with distractor segments flagged as 1
%
% USAGE:
%   bad = flag_distractor_segments(bad, t_dis, tds, dec, MISCorig.Fs);
%
% -------------------------------------------------------------------------

    for n_dis = 1:size(t_dis, 1)
    
        % Convert distractor onset/offset from seconds to audio samples
        n_dis_Y = t_dis(n_dis, :) * Fs;
    
        % Find the closest MEG timepoints to the distractor boundaries
        [~, n_dis_Yds] = min(abs(bsxfun(@minus, tds, n_dis_Y')), [], 2);
    
        % Map to MEG timeline using temporal offset and clamp to valid range
        set_to_bad = dec + n_dis_Yds(1) : dec + n_dis_Yds(2);
        set_to_bad = min(max(set_to_bad, 1), length(bad));
    
        bad(set_to_bad) = 1;
    end
end