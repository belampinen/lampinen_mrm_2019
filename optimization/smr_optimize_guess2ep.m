function ep = smr_optimize_guess2ep(guess, opt)

b               = guess((0 * opt.n_shell + 1):(1 * opt.n_shell));
b_delta         = guess((1 * opt.n_shell + 1):(2 * opt.n_shell));
te              = guess((2 * opt.n_shell + 1):(3 * opt.n_shell));
ndir_scale      = guess((3 * opt.n_shell + 1):end);

% b
b               = smr_optimize_discretize_ep(b,         'b',            opt);

% b_delta
if (isempty(opt.discrete_b_delta))
    b_delta     = smr_optimize_discretize_ep(b_delta,   'b_delta',      opt);
else
    b_delta     = opt.discrete_b_delta(ceil(b_delta));
end

% te
te              = smr_optimize_discretize_ep(te,        'te',           opt);

% ndir
ndir            = opt.discrete_ndir(ceil(ndir_scale));

% EP
ep              = [b(:) b_delta(:) te(:) ndir(:)];

end