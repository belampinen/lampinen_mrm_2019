function s = smr_fit2data_watson(m, xps, my, udirs)
% function s = smr_fit2data_watson(m, xps)
%
% s = s0 * [ f * Ads * At2s + (1 - f) * Adz * At2z]
%
%

% Fiber direction (default)
if (nargin < 3), my = [0 0 1]'; end
if (nargin < 4), uvec_elstat_500; end

%%% Parse input
s0      = m(1);
f       = m(2);
di_s    = m(3);
dd_s    = 1;
di_z    = m(4);
dd_z    = m(5);
od      = m(6);
t2_s    = m(7);
t2_z    = m(8);

%%% Diffusion attenuation
%
% Generate the Watson distribution
kappa = 1 ./ tan((pi / 2) * max(od, eps));
kappa(kappa > 700) = 700;  %
%
w = exp(kappa * (udirs * my).^2);
w = w / sum(w);
%
% For each encoding direction:
% 1) Angle between dir and each point in distribution
alph = xps.u * udirs';
%
% 2) The attenuation correspoding to each point in distribution
f_leg2 = @(x) (3 * x.^2 - 1) / 2;
Ads_full = exp(-di_s * diag(xps.b) * (1 + 2 * dd_s * diag(xps.b_delta) * f_leg2(alph)));
Adz_full = exp(-di_z * diag(xps.b) * (1 + 2 * dd_z * diag(xps.b_delta) * f_leg2(alph)));
%
% 3) The weighted average across distribution
Ads = Ads_full * w;
Adz = Adz_full * w;


%%% Relaxation attenuation
%
At2s = exp(-xps.te/t2_s);
At2z = exp(-xps.te/t2_z);


%%% Signal
%
s = s0 * (f * Ads .* At2s + (1 - f) * Adz .* At2z);

end