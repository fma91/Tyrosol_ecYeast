%% Continue the main-repo G5 dFBA past 93 h (same fit, longer integration)
%
% Calibrates kdScale + uptake boost exactly as tyrosol_dfba (G5, 0–93 h),
% then integrates 0–167 h with those locked parameters.
% Writes dfba_G5_full.* — does NOT overwrite dfba_G5.*.
%
%   addpath('scripts'); tyrosol_dfba_full

extendTo = 167;
tyrosol_dfba;
