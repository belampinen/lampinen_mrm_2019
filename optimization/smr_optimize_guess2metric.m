function metric = smr_optimize_guess2metric(guess, opt)
%
ep              = smr_optimize_guess2ep(guess, opt);
metric          = smr_optimize_ep2metric(ep, opt);
%
if (metric < 0)
    metric = 1e3 * abs(metric);
end
%
end