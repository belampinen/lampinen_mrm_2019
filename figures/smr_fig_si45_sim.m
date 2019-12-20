function smr_fig_si45_sim(output_dir)

if (nargin < 1), output_dir = fullfile(pwd, 'figures', 'output'); end


%%%------------------------------------------------------------------
% PREPARE
%%%------------------------------------------------------------------

do_overwrite = 1;

%
ep_list = {...
    'opt_7_shell', ...
    'opt_6_shell', ...
    'opt_5_shell', ...
    'opt_4_shell', ...    
    'opt_7_shell', ...
    'opt_6_shell', ...
    'opt_5_shell', ...
    'opt_4_shell', ...    
    };
%
mp_name_full_od = 'a_full_od';
mp_name_mid_od  = 'a_mid_od';
mp_list = {...
    mp_name_full_od, ...
    mp_name_full_od, ...
    mp_name_full_od, ...
    mp_name_full_od, ...
    mp_name_mid_od, ...
    mp_name_mid_od, ...
    mp_name_mid_od, ...
    mp_name_mid_od};
c_case_list = 1:8;
n_case = numel(c_case_list);


%%% Simulation
%
opt_optimize            = smr_optimize_opt_derive(smr_optimize_opt);

%
opt_fit                 = smr_opt;
opt_fit.n_rep           = 2;
opt_fit.do_fix_t2       = 0;
%
opt_sim.opt_fit         = opt_fit;
opt_sim.n_realiz        = 200;
opt_sim.do_watson       = 0;
opt_sim.do_gaussian     = 1;




%%%------------------------------------------------------------------
% PERFORM
%%%------------------------------------------------------------------

output_fn = fullfile(output_dir, 'smr_fig_si45_sim.mat');
if (exist(output_fn, 'file') && ~do_overwrite)
    disp(['Exists: ' output_fn])
    return;
end

%%% Simulations
%
sim_mp       = zeros(n_case, opt_sim.n_realiz, 12);
%
for c_c = 1:n_case
    
    c_case = c_case_list(c_c);
    
    ep_name = ep_list{c_case};
    mp_name = mp_list{c_case};
    %
    disp(['ep = ' ep_name ', mp = ' mp_name]);
    disp(' ');
    
    % MP
    mp = smr_mp(mp_name);
        
    % XPS/adjust SNR
    xps = my_ep2xps(smr_ep(ep_name));
    tpf = smr_optimize_xps2tacq(xps, opt_optimize) / opt_optimize.tacq_limit;
    opt_sim.snr = opt_optimize.snr / sqrt(tpf);          
    
     %%% Local minima-robust simulaton
    % Perform twice, using median of first
    %
    for c = 1:2
        if (c == 1)
            opt_sim.opt_fit.init_guess = [];
            sim_mp1 = smr_simest(xps, mp, opt_sim);
        else
            opt_sim.opt_fit.init_guess = median(sim_mp1(:,1:12)) .* [1 1 1e9 1e9 1 1 1 1 1 1 1e3 1e3];
            sim_mp2 = smr_simest(xps, mp, opt_sim);
        end
    end
    %
    sim_mp   (c_case, :, :) = sim_mp2;  
end

% Save
%
save(output_fn, 'sim_mp')
disp(['Wrote: ' output_fn])

end
