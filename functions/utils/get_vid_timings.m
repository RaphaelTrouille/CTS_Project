function [t vid_en_in_SiN t_dis] = get_vid_timings(order,vid)


% img         1   2   3   4   5   6   7   8   9  10  11  12  13  14 
imgtime = {[  0  29  72 100 114 143 160 196 222 262 291 320 354 373 inf],       % 1
           [  0  22  30  45  50  59  78 116 161 202 248 283 316 inf],           % 2
           [  0  24  47  61  89 111 147 189 232 247 264 295 inf]                % 3
           [  0  32  66  99 123 148 159 178 201 231 259 279 inf]                % 4
           [  0  22  63 100 111 161 180 216 267 283 307 331 359 380 400 inf]    % 5 Singing in the file. Cut these periods
           [  0  26  54  96 124 151 203 249 282 319 346 372 inf]                % 6
           [  0  19  50  85 114 145 177 216 245 258 278 296 328 inf]            % 7 Camera bouge parfois, tache sur la camera (toutes les 3 vids).
           [  0  19  56  75  90 109 167 192 199 229 251 294 322 inf]            % 8 
           [  0  13  55  94 141 199 241 282 302 319 344 364 404 inf]            % 9
           [  0  15  28  47  66  80 111 144 165 182 200 218 235 252 285 313 350 372 inf] % 10
           [  0  16  49  81 118 142 154 184 216 267 298 334 inf]                % 11
           [  0  25  46  57  78 103 130 166 197 211 232 254 268 inf]};          % 12
% onsets and offsets of the videos [s]
onset = [13 389     % 1
         10 356     % 2
         11 387     % 3
         24 323     % 4
          9 442     % 5
          9 409     % 6  Caution. Video ends at 290 s
          2 363     % 7
          3 333     % 8
          2 439     % 9
          4 391     % 10
          7 386     % 11
         11 279];   % 12
% timing to keep in the video but not used in the analysis (singing)
disregardtime = {[],       % 1
                 [],       % 2
                 [],       % 3
                 [],       % 4
                 [131 148; 187 203; 238 256; 425 442],       % singing
                 [],       % 6
                 [238 245],       % 7
                 [],       % 8
                 [],       % 9
                 [],       % 10
                 [],       % 11
                 []};      % 12
%                 1 2 3 4 5 6 7 8 9 10    11: first 5 seconds
cond_vid =       [0 0 0 1 0 1 0 1 0 1 1];
cond_energetic = [0 0 0 0 1 1 0 0 1 1 0];
cond_info =      [0 0 0 0 0 0 1 1 1 1 0];
cond_SiN =       [0 0 1 1 1 1 1 1 1 1 0];




all_pick_vid = [1 4 9 12
                2 6 7 11
                3 5 8 10];

all_orders =[  3     8     9     6     1     4    10     7     5     2
               1     9     6     3     8     2     4    10     7     5
               9     3     1     8     6    10     2     5     4     7
               6     1     8     9     3     5     7     4     2    10];
all_orders([2 3],:) = all_orders([2 3],[6:10 1:5]);

switch order
    case 1
    case 2
        all_orders = all_orders([2 1 4 3],:);
    case 3
        all_orders = all_orders([3 4 1 2],:);
    case 4
        all_orders = all_orders([4 3 2 1],:);
end

% first 5 seconds: no noise, video
all_orders = cat(2,ones(4,1)*11,all_orders);


% Video 6 was does not contain the whole narrative. Force the 3 last
% conditions to be non-video ones.
if vid == 6
    this_order = all_orders(2,7:11);
    [dummy,new_order] = sort(cond_vid(this_order),'descend');
    all_orders(2,7:11) = this_order(new_order);
end


n_vid = ceil(vid/3);
all_orders = all_orders(n_vid,:);




% get the timing of the conditions and images, taking into
% account the time to disregard (but keep in the video)
t_dis = disregardtime{vid}; L_dis = sum(diff(t_dis,[],2));
t = [onset(vid,1) linspace(onset(vid,1)+5,onset(vid,2)-L_dis,11)];
for k_dis = 1:size(t_dis,1)
    ind = find(t < t_dis(k_dis,1),1,'last');
    t(ind+1:end) = t(ind+1:end) + diff(t_dis(k_dis,:));
end
t_img = imgtime{vid};

all_orders = all_orders(2:end);
t_dis = t_dis-t(1);
t = t(2:end)-t(1);


vid_en_in_SiN = [cond_vid ; cond_energetic ; cond_info ; cond_SiN];
vid_en_in_SiN = vid_en_in_SiN(:,all_orders);

