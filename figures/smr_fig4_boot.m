function [mp_fit_vec, s_vec, s_fit_vec] = smr_fig4_boot(xps, mp, opt)

if (~isfield(opt, 'do_watson'))  , opt.do_watson   = 0; end
if (~isfield(opt, 'do_gaussian')), opt.do_gaussian = 0; end

%%% PREPARE
%
sigma_noise     = 1 / opt.snr;  % S0 normalized to 1
%
if (opt.do_watson)
    s               = smr_fit2data_watson(mp, xps, [0 0 1]', uvec_elstat_500);
else
    s               = dtd_fit2data(mp, xps);
end
% 

%%% PERFORM
%
mp_fit_vec      = zeros(opt.n_boot, 13);
s_vec           = zeros(opt.n_boot, xps.n);
s_fit_vec       = s_vec;
%
for c_boot = 1:opt.n_boot

    disp(['c_boot = ' num2str(c_boot)])
    
    % Noise
    if (opt.do_gaussian)
        s_noised = s + randn([xps.n 1]) * sigma_noise;
    else
        s_noised = sqrt( (s + randn([xps.n 1]) * sigma_noise).^2 + (0 + randn([xps.n 1]) * sigma_noise).^2 );
    end
    
    % Fit
    mp_fit   = smr_data2fit(s_noised, xps, opt.opt_fit);
    s_fit    = smr_fit2data(mp_fit, xps);
    
    % Store
    s_vec      (c_boot, :)  = s_noised;
    s_fit_vec  (c_boot, :)  = s_fit;
    mp_fit_vec (c_boot, :)  = mp_fit; 
end
end