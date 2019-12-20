function smr_fig_si1_sim(output_dir)

if (nargin < 1), output_dir = fullfile(pwd, 'figures', 'output'); end

%%%------------------------------------------------------------------
% PREPARE
%%%------------------------------------------------------------------

do_overwrite            = 1;

ep_list = {...
    'opt_lte', ...                 % 1
    'opt_dde', ...                 % 2
    'opt_ste', ...                 % 3
    ...
    'opt_4_shell', ...             % 4
    'opt_5_shell', ...             % 5
    'opt_6_shell', ...             % 6
    'opt_7_shell', ...             % 7
    'opt_8_shell', ...             % 8
    'opt_9_shell', ...             % 9
    'opt_10_shell', ...            % 10
    };
%
c_ep_process = 1:10;
n_ep = numel(c_ep_process);


%%% MP
%
%
mp_list = {'a', 'b', 'c'};
%
c_mp_process = 1:3;
n_mp = numel(c_mp_process);
%


%%% Simulation
%
opt_optimize            = smr_optimize_opt_derive(smr_optimize_opt);
opt_optimize.par_sign_diff(5) = 0.1;
snr                     = opt_optimize.snr;

%
opt_fit                 = smr_opt;
opt_fit.n_rep           = 2;
opt_fit.do_fix_t2       = 0;
%
opt_sim.n_realiz        = 100;
opt_sim.snr             = snr;
opt_sim.opt_fit         = opt_fit;
opt_sim.do_watson       = 0;



%%%------------------------------------------------------------------
% PERFORM
%%%------------------------------------------------------------------


output_fn = fullfile(output_dir, 'smr_fig_si1_sim.mat');
if (exist(output_fn, 'file') && ~do_overwrite)
    disp(['Exists: ' output_fn])
    return;
end
%
vw           = zeros(n_ep, 1);
v_vec        = zeros(n_ep, n_mp, 12);
%
for c_ep = 1:n_ep

    disp(['ep = ' ep_list{c_ep_process(c_ep)}])
    %
    ep  = smr_ep(ep_list{c_ep_process(c_ep)});
    xps = my_ep2xps(ep);   
    t_acq = smr_optimize_xps2tacq(xps, opt_optimize);
        
    % Estimate variances
    %
    for c_mp = 1:n_mp
        
        disp(['mp = ' mp_list{c_mp_process(c_mp)}]);
        disp(' ');
                
        mp = smr_mp(mp_list{c_mp_process(c_mp)});
        
        %%% Local minima-robust simulaton
        % Perform twice, using median of first
        %
        for c = 1:2
            if (c == 1)
                opt_sim.opt_fit.init_guess = [];
                mp1 = smr_simest(xps, mp, opt_sim);
            else
                opt_sim.opt_fit.init_guess = median(mp1(:,1:12)) .* [1 1 1e9 1e9 1 1 1 1 1 1 1e3 1e3];
                mp2 = smr_simest(xps, mp, opt_sim);
            end
        end
        %
        v_vec(c_ep,c_mp,:) = var(mp2);    
    end
    
   % Convert variances to vw
   v_tmp = squeeze(v_vec(c_ep, :, :));
   %
   w        = (1 ./ opt_optimize.par_sign_diff.^2)';
   tmp      = mean(v_tmp * w);
   vw(c_ep) = tmp * (t_acq / opt_optimize.tacq_limit);
end

save(output_fn, 'vw')
disp(['Wrote: ' output_fn]);
end
