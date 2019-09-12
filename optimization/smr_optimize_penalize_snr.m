function snr_pf = smr_optimize_penalize_snr(snr, snr_limit)
%
do_penalize         = snr < snr_limit;
snr_pf              = (50 * (snr - snr_limit)).^2 .* do_penalize + 1;
end