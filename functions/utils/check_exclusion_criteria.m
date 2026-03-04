function [should_skip, reason] = check_exclusion_criteria(sub_name, bad_ica, bad_diagnosis)
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