function smr_fig_3_sim(output_dir)
%
% Sweeps a locking of fs from 0 to 1 using ep
% and adding noise to ground truth mp
%
if (nargin < 1), output_dir = fullfile(pwd, 'figures', 'output'); end


%%%------------------------------------------------------------------
% PREPARE
%%%------------------------------------------------------------------

do_overwrite = 1;

%%% EP and MP
ep_list = {...
    'opt_10_shell', ...
    'in_vivo', ...
    'opt_lte', ...
    'lampinen_hbm_2019'};
n_ep = numel(ep_list);
%
mp_list = {'a', 'b', 'c'};
n_mp    = numel(mp_list);

%%% NRV iteration
%
opt_optimize        = smr_optimize_opt_derive(smr_optimize_opt);
%
opt_fit             = smr_opt;
opt_fit.n_rep       = 2;
opt_fit.do_fix_t2   = 0;
%
n_iter              = 10;
step                = 0.0125;
fs_fix              = (step/2):step:(1-step/2);
n_fix               = numel(fs_fix);
%

%%%------------------------------------------------------------------
% PERFORM
%%%------------------------------------------------------------------

%
for c_mp = 1:n_mp
    %
    mp = smr_mp(mp_list{c_mp});
    %
    for c_ep = 1:n_ep                
        
        output_fn = fullfile(output_dir, ['smr_fig_3_sim_' mp_list{c_mp} '_' ep_list{c_ep} '.mat']);
        if (exist(output_fn, 'file') && ~do_overwrite)
            disp(['Exists: ' output_fn]);
            continue;
        end
        %
        disp(['mp = ' mp_list{c_mp}, ', ep = ' ep_list{c_ep}]);
        disp(' ');
        %
        
        % Adjust SNR
        % 
        ep = smr_ep(ep_list{c_ep});
        xps = my_ep2xps(ep);
        tacq_factor = smr_optimize_xps2tacq(xps, opt_optimize) / opt_optimize.tacq_limit;
        my_snr = opt_optimize.snr / sqrt(tacq_factor);
        
        % Create signal
        %
        s           = smr_fit2data(mp, xps);
        sigma_noise = 1 / my_snr; % S0 normalized to 1

        % Iterate
        %
        ssr_vals = zeros(n_iter, n_fix);
        mp_vals  = zeros(n_iter, n_fix, 12);
        %
        for c_iter = 1:n_iter
            
            disp(['c_iter = ' num2str(c_iter)]);
            
            
            % Noise
            s_noised = s + randn([xps.n 1]) * sigma_noise;
            
            % Fit with different fixed values for fs
            for c_fix = 1:n_fix
                
                % Fix fs
                opt_fit.f_min    = fs_fix(c_fix) - step / 2;
                opt_fit.f_max    = fs_fix(c_fix) + step / 2;
                
                % Fit
                mp_fit  = smr_data2fit(s_noised, xps, opt_fit);
                
                % Evaluate
                s_fit   = smr_fit2data(mp_fit, xps);
                ssr     = sum((s_noised - s_fit).^2);
                
                % Store
                ssr_vals(c_iter, c_fix)    = ssr;
                mp_vals (c_iter, c_fix, :) = mp_fit(1:12);
            end
        end  
        % Save
        %
        save(output_fn, 'ssr_vals', 'mp_vals', 'xps')
        disp(['Wrote: ' output_fn]);
    end
end
end