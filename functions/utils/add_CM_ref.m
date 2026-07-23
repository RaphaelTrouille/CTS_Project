function CM = add_CM_ref(CM, info, chan, filt, norm_val)
% ADD_CM_REF  Add a reference signal to a CM (Cortical Mapping) structure.
%
% DESCRIPTION:
%   Appends a new reference signal entry to the CM.ref array. Each reference
%   is described by a label, a signal vector, an optional filter specification,
%   and an optional normalization flag. All references in a CM structure must
%   have the same number of samples — a consistency check is enforced when
%   adding subsequent references.
%
% INPUTS:
%   CM        - CM structure with field .ref (may be empty for the first entry)
%   info      - Label/description string for the reference (e.g. 'global sound')
%   chan      - Signal vector for the reference channel (1 x N)
%   filt      - (Optional) Filter specification. Default: [] (no filtering)
%   norm_val  - (Optional) Normalization flag or value. Default: 0 (no normalization)
%
% OUTPUTS:
%   CM  - Updated CM structure with the new reference appended to CM.ref
%         Each entry in CM.ref contains fields:
%           .info : label string
%           .chan : signal vector
%           .filt : filter specification
%           .norm : normalization value
%
% USAGE:
%   CM = add_CM_ref(CM, 'global sound', env_signal);
%   CM = add_CM_ref(CM, 'Mouth Surface', Vmouth_FS, [], 1);
%
% -------------------------------------------------------------------------    

    % Set default values for optional arguments
    if nargin < 4, filt = []; end
    if nargin < 5, norm_val = 0; end
    
    if isempty(CM.ref)
        % First reference: no consistency check needed
        Nsig = 1;
    else
        % Subsequent references: enforce same length as existing ones
        if length(chan) ~= length(CM.ref(1).chan)
            error('Error in add_CM_ref: Signal "%s" has not the same length as previous references', info);
        end
        Nsig = length(CM.ref) + 1;    
    end

    % Append new reference entry
    CM.ref(Nsig).info = info;
    CM.ref(Nsig).chan = chan;
    CM.ref(Nsig).filt = filt;
    CM.ref(Nsig).norm = norm_val;

    fprintf('Ref #%d [%s] added successfully.\n', Nsig, info);
end