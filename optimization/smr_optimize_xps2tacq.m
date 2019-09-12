function [tacq, tr] = smr_optimize_xps2tacq(xps, opt)

g       = smr_optimize_xps2g(xps, opt);
%
tseq    = opt.tau_fatsat + max(xps.te) + (opt.tau_epi / 2);
%
tr      = max(opt.duty_cycle * max(g), tseq * opt.n_slice);
tacq    = xps.n * tr;

end