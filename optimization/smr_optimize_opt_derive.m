function opt = smr_optimize_opt_derive(opt)

%%% Bounds
%
b_bounds        = repmat([opt.b_min          opt.b_max],          [opt.n_shell 1]);
%
te_bounds       = repmat([opt.te_min         opt.te_max],         [opt.n_shell 1]);
%
if (isempty(opt.discrete_b_delta))
    b_delta_bounds  = repmat([opt.b_delta_min    opt.b_delta_max],[opt.n_shell 1]);
else
    b_delta_bounds =  repmat([0 numel(opt.discrete_b_delta)], [opt.n_shell 1]);
end
%
ndir_bounds     = repmat([eps numel(opt.discrete_ndir)], [opt.n_shell 1]);
%
opt.bounds = [b_bounds; b_delta_bounds; te_bounds; ndir_bounds];


%%% Sequence timing (unbalanced gradients)
%
opt.tau_epi_prim            = max(opt.tau_epi * (opt.pf_factor - 0.5), 0);
te_minimal_left             = (0.5 * opt.tau_rf90); % Half RF pulse 
te_minimal_right            = opt.tau_epi_prim;    
%
opt.te_minimal              = te_minimal_left + te_minimal_right + opt.tau_rf180;
%
% Time not available for diffusion encoding
opt.te_nonenc = opt.te_minimal - opt.tau_rf180;


%%% SNR
%
opt.snr       = opt.snr_ref * prod(opt.voxel_size) / prod(opt.snr_ref_voxel_size);

end