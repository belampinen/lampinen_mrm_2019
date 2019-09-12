function [ssr_vals, mp_vals] = smr_fig3_iterate(xps, mp, opt)
%
% Sweeps a locking of fs from 0 to 1 using ep
% and adding noise to ground truth mp
%
if (~isfield(opt, 'do_watson'))  , opt.do_watson   = 0; end
if (~isfield(opt, 'do_gaussian')), opt.do_gaussian = 0; end


%%% PREPARE
%
if (opt.do_watson)
    s           = smr_fit2data(mp, xps, [0 0 1]', uvec_elstat_500);
else
    s           = smr_fit2data(mp, xps);
end
%
sigma_noise     = 1 / opt.snr; % S0 normalized to 1
%


%%% PERFORM
%
if (isfield(opt, 'fs'))
    par_fix = opt.fs;
elseif (isfield(opt, 'dd_z'))
    par_fix = opt.dd_z;
else
    step    = 0.025;
    par_fix = (step/2):step:(1-step/2);
end
step  = par_fix(2) - par_fix(1);
n_fix = numel(par_fix);

% Iterate
ssr_vals = zeros(opt.n_iter, n_fix);
mp_vals  = zeros(opt.n_iter, n_fix, 13);
%
for c_iter = 1:opt.n_iter
    
    disp(['c_iter = ' num2str(c_iter)]);
    
    % Noise
    if (opt.do_gaussian)
        s_noised = s + randn([xps.n 1]) * sigma_noise;
    else
        s_noised = sqrt( (s + randn([xps.n 1]) * sigma_noise).^2 + (0 + randn([xps.n 1]) * sigma_noise).^2 );
    end
    
    % Introduce outliers
    if (isfield(opt, 'ol_freq'))
        s_noised = sqrt( s_noised.^2 + 5 * sigma_noise * randn(size(s)).^2 .* (rand(size(s)) < opt.ol_freq) );
    end
    
    for c_fix = 1:n_fix
        
        % Fix fs
        opt.opt_fit.f_min    = par_fix(c_fix) - step / 2;
        opt.opt_fit.f_max    = par_fix(c_fix) + step / 2;
        
        % Fit
        mp_fit  = smr_data2fit(s_noised, xps, opt.opt_fit);
        
        % Evaluate
        s_fit   = smr_fit2data(mp_fit, xps);
        ssr     = sum((s_noised - s_fit).^2);
        
        % Store
        ssr_vals(c_iter, c_fix)    = ssr;
        mp_vals (c_iter, c_fix, :) = mp_fit;
        
    end
end
end