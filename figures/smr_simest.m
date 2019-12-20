function [mp_fit_vec, s_vec, s_fit_vec] = smr_simest(xps, mp, opt_sim)

if (~isfield(opt_sim, 'do_watson'))  , opt_sim.do_watson   = 0; end
if (~isfield(opt_sim, 'do_gaussian')), opt_sim.do_gaussian = 0; end

%%% PREPARE
%
sigma_noise     = 1 / opt_sim.snr;  % S0 normalized to 1
%
if (opt_sim.do_watson)
    s               = smr_fit2data_watson(mp, xps, [0 0 1]', uvec_elstat_500);
else
    s               = smr_fit2data(mp, xps);
end
% 

%%% PERFORM
%
mp_fit_vec      = zeros(opt_sim.n_realiz, 12);
s_vec           = zeros(opt_sim.n_realiz, xps.n);
s_fit_vec       = s_vec;
%
for c_realiz = 1:opt_sim.n_realiz

    disp(['c_realiz = ' num2str(c_realiz)])
    
    % Noise
    if (opt_sim.do_gaussian)
        s_noised = s + randn([xps.n 1]) * sigma_noise;
    else
        s_noised = sqrt( (s + randn([xps.n 1]) * sigma_noise).^2 + (0 + randn([xps.n 1]) * sigma_noise).^2 );
    end
    
    % Fit
    mp_fit   = smr_data2fit(s_noised, xps, opt_sim.opt_fit);    
    s_fit    = smr_fit2data(mp_fit, xps);
    
    % Store
    s_vec      (c_realiz, :)  = s_noised;
    s_fit_vec  (c_realiz, :)  = s_fit;
    mp_fit_vec (c_realiz, :)  = mp_fit(1:12); 
end
end