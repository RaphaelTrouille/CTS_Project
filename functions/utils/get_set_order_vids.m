function [set order vids] = get_set_order_vids(subfold)

files = dir(fullfile(subfold,'*.fif'));
Nfiles = length(files);
set = ones(1,Nfiles)*NaN;
order = ones(1,Nfiles)*NaN;
vids = ones(1,Nfiles)*NaN;
for n_file = 1:length(files)
    this_file = files(n_file).name;
    if length(this_file) > 26
        set(n_file) = str2double(this_file(13));
        order(n_file) = str2double(this_file(20));
        vids(n_file) = str2double(this_file(25:26));
        if isnan(vids(n_file))
            vids(n_file) = str2double(this_file(25));
        end
    end
end
set = unique(set(~isnan(set)));
order = unique(order(~isnan(order)));
vids = unique(vids(~isnan(vids)));

% Check if set order were available
if isempty(vids)
    set = ones(1,Nfiles)*NaN;
    order = ones(1,Nfiles)*NaN;
    vids = ones(1,Nfiles)*NaN;
    for n_file = 1:length(files)
        this_file = files(n_file).name
        if strcmp(this_file(10:12),'vid')
            vids(n_file) = str2double(this_file(13:14));
            if isnan(vids(n_file))
                vids(n_file) = str2double(this_file(13));
            end
            if any(not(vids(n_file)-[1 4 9 12]))
                set(n_file) = 1;
            elseif  any(not(vids(n_file)-[2 6 7 11]))
                set(n_file) = 2;
            else
                set(n_file) = 3;
            end
            order(n_file) = 1;
        end
    end
end
set = unique(set(~isnan(set)));
order = unique(order(~isnan(order)));
vids = unique(vids(~isnan(vids)));

if length(set) ~= 1
    error('number of sets not equal to 1')
end
if length(order) ~= 1
    error('number of orders not equal to 1')
end

if length(vids) ~= 4
    error('number of vids not equal to 4')
end