function out = iif(cond, valTrue, valFalse)
% IIF  Inline if — return one of two values based on a condition.
%
% DESCRIPTION:
%   Functional equivalent of a ternary operator. Returns valTrue if cond
%   is true, valFalse otherwise. Useful for concise inline assignments
%   without writing a full if/else block.
%
% INPUTS:
%   cond      - Scalar logical condition to evaluate
%   valTrue   - Value returned if cond is true
%   valFalse  - Value returned if cond is false
%
% OUTPUTS:
%   out  - Either valTrue or valFalse depending on cond
%
% USAGE:
%   label = iif(score > 0.5, 'good', 'bad');
%   n     = iif(isempty(x), 0, length(x));
%
% NOTE:
%   Both valTrue and valFalse are evaluated before the condition is checked
%   (eager evaluation). Avoid using this with expressions that have side
%   effects or are computationally expensive.
%
% -------------------------------------------------------------------------
 
    if cond, out = valTrue; else, out = valFalse; end
end