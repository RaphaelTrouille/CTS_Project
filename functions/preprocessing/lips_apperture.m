function Vmouth_FS = lips_apperture(vid_dir, vid_n, attended_speech)

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