function opt = smr_opt(opt)
% function opt = smr_opt(opt)
%
% Makes sure that all needed fields in the options structure are present
%
% Adapted for the MDM-structure

opt.present = 1;

opt = my_ensure_field(opt, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off','MaxFunEvals',1e4));

opt = my_ensure_field(opt, 'fig_maps', ...
    {'s0', ...
    'fs', 'di_s', 'di_z', 'dd_z', ...
    't2_s', 't2_z', ...
    'p2', 'da_s', 'da_z', 'dr_z', 't2', 'ssr1', 'ssr2', 'ssr3'});

opt = my_ensure_field(opt, 'do_fix_t2', 0);
opt = my_ensure_field(opt, 'do_fix_da', 0);
opt = my_ensure_field(opt, 'do_tortuosity', 0);

opt = my_ensure_field(opt, 'fig_prefix', 'smr');

opt = my_ensure_field(opt, 'n_rep', 1);

opt = my_ensure_field(opt, 'init_guess', []);

