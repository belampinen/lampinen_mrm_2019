function [of, export] = smr_optimize_ep2metric(ep, opt)

%%% MP
mp      = smr_optimize_mp(opt.mp_name);
n_mp    = size(mp, 1);

%%% XPS
xps     = my_ep2xps(ep);
t_acq   = smr_optimize_xps2tacq(xps, opt);

%%% Vw
vw    = smr_optimize_ep2vw(ep, opt);

%%% OF (apply penalties to Vw)
%
% Noise floor (consider each mp)
snr_pf = 1;
if (~isempty(opt.snr_min))
    snr_min = zeros(n_mp, 1);
    for c_mp = 1:n_mp
        snr_min(c_mp) = min(opt.snr * my_pa(smr_fit2data(mp(c_mp,:), xps), xps));
    end
    snr_pf = smr_optimize_penalize_snr(min(snr_min), opt.snr_min);
end
%
% Gradient
g       = smr_optimize_xps2g(xps, opt);
g_pf    = smr_optimize_penalize_g(max(g), opt.g_max);
%
%
of  = vw * snr_pf * g_pf;

%%% Export
export.vw       = vw;
export.t_acq    = t_acq;
export.g_max    = max(g);
export.g_pf     = g_pf;
export.snr_min  = min(snr_min);
export.snr_pf   = snr_pf;

end


