function opt = smr_opt(opt)
% function opt = smr_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.present = 1;

opt = my_ensure_field(opt, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off','MaxFunEvals',1e4));

opt = my_ensure_field(opt, 'fig_maps', ...
    {'s0', ...
    'fs', 'di_s', 'di_z', 'dd_z', ...
    't2_s', 't2_z', ...
    'msr', ...
    'p2', 'da_s', 'da_z', 'dr_z', 't2'});

opt = my_ensure_field(opt, 'fig_prefix', 'smr');

opt = my_ensure_field(opt, 'n_rep', 1);

opt = my_ensure_field(opt, 'init_guess', []);

