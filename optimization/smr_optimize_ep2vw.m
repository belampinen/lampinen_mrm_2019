function vw = smr_optimize_ep2vw(ep, opt)

if (nargin < 2)
    opt = smr_optimize_opt_derive(smr_optimize_opt);
end

%%% XPS
xps     = my_ep2xps(ep);
t_acq   = smr_optimize_xps2tacq(xps, opt);

%%% MP
mp      = smr_optimize_mp(opt.mp_name);
n_mp    = size(mp, 1);

%%% CRLB
%
crlb  = zeros(size(mp));
for c_mp = 1:n_mp
    crlb(c_mp,:) = smr_crlb(xps, mp(c_mp,:), opt.snr);
end

%%% Weighting and averaging
%
% Weight CRBL by inverse square significant difference
w       = (1 ./ opt.par_sign_diff.^2)';
vw      = mean(crlb * w);

%%% T_acq factor
% Scale Vw as if tacq = tacq_limit, assuming SNR ~ sqrt(tacq)
%
vw   = vw * (t_acq / opt.tacq_limit);

end