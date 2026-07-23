function Vmouth_FS = lips_apperture(vid_dir, vid_n, attended_speech)
% LIPS_APPERTURE  Estimate lip aperture from facial landmark coordinates
%                 and resample it to the audio sampling rate.
%
% DESCRIPTION:
%   This function reads facial landmark coordinates extracted from a video
%   file, computes the mouth aperture signal as the signed area difference
%   between the upper and lower lip contours (using the shoelace formula),
%   applies a high-pass cosine filter to remove slow drift, and finally
%   resamples the signal from the video frame rate (Fm = 50 Hz) to the
%   audio sampling rate (FS) using Gaussian kernel interpolation.
%
% INPUTS:
%   vid_dir         - Path to the directory containing landmark coordinate files
%   vid_n           - Video number (integer), used to build the filename
%   attended_speech - Structure with fields:
%                       .signal : audio signal vector (used for sizing output)
%                       .Fs     : audio sampling rate in Hz
%
% OUTPUTS:
%   Vmouth_FS  - Lip aperture signal resampled at the audio sampling rate (FS),
%                same length as attended_speech.signal
%
% FILE FORMAT:
%   The function expects a text file named:
%     'Vid<vid_n>_SpeechTrack_face_coordinates.txt'
%   in vid_dir, containing a matrix of facial landmark coordinates.
%   Each landmark k has x and y stored at columns (k*8-7) and (k*8-6).
%
% METHOD:
%   - Upper lip contour: landmarks [4 6 9 1 10 5 3]
%   - Lower lip contour: landmarks [4 8 2 7 3]
%   - Aperture = signed area(upper) - signed area(lower), via trapezoidal rule
%   - High-pass filtered at 0.1 Hz to remove head-motion drift
%   - Resampled to FS via Gaussian kernel (sigma = 0.6 * FS/Fm)
%
% DEPENDENCIES:
%   - cosine_filter.m (custom function, must be on MATLAB path)
%
% USAGE:
%   Vmouth_FS = lips_apperture(vid_dir, 3, attended_speech);
%
% -------------------------------------------------------------------------

file_path = fullfile(vid_dir, ['Vid' num2str(vid_n) '_SpeechTrack_face_coordinates.txt']);

FS = attended_speech.Fs;

if ~exist(file_path, 'file')
    error('File does not exist: %s', file_path);
end

pnts = readmatrix(file_path);
Fm = 50;

% Upper lip
keep_u = [4 6 9 1 10 5 3];
x_ulip = pnts(:, keep_u * 8 - 7);
y_ulip = pnts(:, keep_u * 8 - 6);

% Lower Lip
keep_l = [4 8 2 7 3];
x_llip = pnts(:, keep_l * 8 - 7);
y_llip = pnts(:, keep_l * 8 - 6);

Vup = sum((y_ulip(:, 1:end-1) + y_ulip(:, 2:end)) .* diff(x_ulip, [], 2), 2) / 2;
Vlow = sum((y_llip(:, 1:end-1) + y_llip(:, 2:end)) .* diff(x_llip, [], 2), 2) / 2;
Vmouth = Vup' - Vlow';

% Filtering
F_h_x = cosine_filter(length(Vmouth), {'high'}, 0.1/Fm*2, 0.1/Fm*2);
Vmouth = real(...
            ifft(...
                fft(Vmouth) .* F_h_x));
Vmouth_FS = zeros(size(attended_speech.signal));
Vspread = zeros(size(attended_speech.signal));

t_frames = round(Fm/2 : FS/Fm:length(attended_speech.signal));

n_samples = min(length(t_frames), length(Vmouth));
idx_target = t_frames(1:n_samples);

Vmouth_FS(idx_target) = Vmouth(1:n_samples);
Vspread(idx_target) = 1;

sigma = 0.6 * FS / Fm;
x_kern = -3*FS/Fm : 3*FS/Fm;
kern = normpdf(x_kern, 0, sigma);

Vmouth_FS = conv(Vmouth_FS, kern, 'same');
Vspread   = conv(Vspread, kern, 'same');

Vmouth_FS = Vmouth_FS ./ Vspread;

Vmouth_FS(isnan(Vmouth_FS)) = 0;
end