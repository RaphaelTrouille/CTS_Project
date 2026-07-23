function [should_skip, reason] = check_exclusion_criteria(sub_name, bad_ica, bad_diagnosis)
% CHECK_EXCLUSION_CRITERIA  Determine whether a subject should be excluded
%                            from analysis based on predefined exclusion lists.
%
% DESCRIPTION:
%   Checks whether a given subject appears in any exclusion list and returns
%   a flag and reason accordingly. Currently two exclusion criteria are
%   supported: incomplete ICA correction and diagnostic/IQ concerns.
%
% INPUTS:
%   sub_name       - Subject identifier string (e.g. 'meg01')
%   bad_ica        - Cell array of subject IDs with incomplete ICA correction
%   bad_diagnosis  - Cell array of subject IDs with diagnostic or IQ concerns
%
% OUTPUTS:
%   should_skip  - Boolean; true if the subject meets any exclusion criterion
%   reason       - String describing the exclusion reason, or '' if none
%
% USAGE:
%   [skip, reason] = check_exclusion_criteria('meg04', bad_ica, bad_diagnosis);
%   if skip, fprintf('Skipping %s: %s\n', sub_name, reason); end
%
% -------------------------------------------------------------------------    
    
    should_skip = false;
    reason = '';

    if any(contains(bad_ica, sub_name))
        should_skip = true;
        reason = 'ICA not corrected yet';
    elseif any(contains(bad_diagnosis, sub_name))
        should_skip = true;
        reason = 'Diagnosis/IQ doubts';
    end
end