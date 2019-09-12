function [vw, t_acq] = smr_optimize_vw(xps, mp, opt)

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
t_acq = smr_optimize_xps2tacq(xps, opt);
vw   = vw * (t_acq / opt.tacq_limit);


end