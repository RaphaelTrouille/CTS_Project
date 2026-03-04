function CM = add_CM_ref(CM, info, chan, filt, norm_val)
    if nargin < 4, filt = []; end
    if nargin < 5, norm_val = 0; end

    if isempty(CM.ref)
        Nsig = 1;
    else
        if length(chan) ~= length(CM.ref(1).chan)
            error('Error in add_CM_ref: Signal "%s" has not the same lenth as previous references', info);
        end
        Nsig = length(CM.ref) + 1;    
    end

    CM.ref(Nsig).info = info;

    CM.ref(Nsig).chan = chan;
    CM.ref(Nsig).filt = filt;
    CM.ref(Nsig).norm = norm_val;

    fprintf('Ref #%d [%s] added successfully.\n', Nsig, info);
end