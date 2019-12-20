function smr_fig_4_sim(output_dir)

if (nargin < 1), output_dir = fullfile(pwd, 'figures', 'output'); end

%%%------------------------------------------------------------------
% PREPARE
%%%------------------------------------------------------------------

do_overwrite = 1;
%
%
ep_name = 'in_vivo';
%
xps     = my_ep2xps(smr_ep(ep_name));
%
mp_list = {'zero_od', 'mid_od', 'full_od', 'full_od_plus'};
%
n_mp    = numel(mp_list);
%

%
opt_optimize            = smr_optimize_opt_derive(smr_optimize_opt);
snr                     = opt_optimize.snr;
%
opt_fit                 = smr_opt;
opt_fit.n_rep           = 2;
opt_fit.do_fix_t2       = 0;
%
opt_sim.n_realiz        = 200;
opt_sim.snr             = snr;
opt_sim.opt_fit         = opt_fit;
opt_sim.do_watson       = 1;



%%%------------------------------------------------------------------
% PERFORM
%%%------------------------------------------------------------------


output_fn = fullfile(output_dir, 'smr_fig_4_sim.mat');
if (exist(output_fn, 'file') && ~do_overwrite)
    disp(['Exists: ' output_fn])
    return;
end


%%% Simulations
%
sim_mp       = zeros(n_mp, opt_sim.n_realiz, 12);
sim_s        = zeros(n_mp, opt_sim.n_realiz, xps.n);
sim_s_fit    = sim_s;
%
for c_mp = 1:n_mp
    
    disp(['mp = ' mp_list{c_mp}]);
    disp(' ');
    %
    mpw = smr_mp_watson(mp_list{c_mp});
    
    
    %%% Local minima-robust simulaton
    % Perform twice, using median of first
    %
    for c = 1:2
        if (c == 1)
            opt_sim.opt_fit.init_guess = [];
            [sim_mp1, ~, ~]               = smr_simest(xps, mpw, opt_sim);
        else
            opt_sim.opt_fit.init_guess = median(sim_mp1(:,1:12)) .* [1 1 1e9 1e9 1 1 1 1 1 1 1e3 1e3];
            [sim_mp2, sim_s2, sim_s_fit2] = smr_simest(xps, mpw, opt_sim);
        end
    end
    %
    sim_mp   (c_mp, :, :) = sim_mp2;
    sim_s    (c_mp, :, :) = sim_s2;
    sim_s_fit(c_mp, :, :) = sim_s_fit2;
end

% Save
%
save(output_fn, ...
    'sim_mp', 'sim_s', 'sim_s_fit', 'xps')
disp(['Wrote: ' output_fn]);

end