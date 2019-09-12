function smr_fig4(output_dir)

if (nargin < 1), output_dir = fullfile(pwd, 'figures', 'output'); end

do_overwrite            = 0;
output_name             = 'smr_fig4.tiff';


%%%------------------------------------------------------------------
% PREPARE
%%%------------------------------------------------------------------


%%% Simulation
%
%
ep_name     = 'in_vivo';
xps         = my_ep2xps(smr_ep(ep_name));
%
mp_list = {'zero_od', 'mid_od', 'full_od', 'full_od_plus'};
n_mp    = numel(mp_list);
%
%
opt_optimize            = smr_optimize_opt_derive(smr_optimize_opt);
snr                     = opt_optimize.snr;
%
%
opt_fit                 = smr_opt;
opt_fit.n_rep           = 2;
%
opt_boot.n_boot         = 200;
opt_boot.snr            = snr;
opt_boot.opt_fit        = opt_fit;
opt_boot.do_watson      = 1;


%%% Plotting
%
xticklbl    = {'Zero OD', 'Mid OD', 'Full OD', 'Full OD+'};
%
c_par_list   = [2 3 4 5 11 12 14]; % fs, dis, diz, ddz, T2s, T2z, p2
c_par_list_w = c_par_list;
c_par_list_w(c_par_list_w > 5) = c_par_list_w(c_par_list_w > 5) - 4;
n_par        = numel(c_par_list);
%
n_row          = 2;
n_col          = 4;
%
par_list = {...
    'S0', ...
    '\itf\rm\bf_{S}', ...
    '\itD\rm\bf_{I;S} [\mum^2/ms]', ...
    '\itD\rm\bf_{I;Z} [\mum^2/ms]', ...
    '\itD\rm\bf_{\Delta;Z}', ...
    [], ...
    [], ...
    [], ...
    [], ...
    [], ...
    'T_{2;S} [ms]', ...
    'T_{2;Z} [ms]', ...
    [], ...
    '\itp\rm\bf_2'};
%
par_sc = [1 1 1e9 1e9 1 1 1 1 1 1 1e3 1e3 1 1];
%
par_lims = {...
    [], ...
    [0 1], ...
    [0 1.5], ...
    [0 1.5], ...
    [-0.5 1], ...
    [], ...
    [], ...
    [], ...
    [], ...
    [], ...
    [50 100], ...
    [50 100], ...
    [0 300], ... % msr
    [0 1]
    };
par_tick_lbl = {...
    [], ...
    {'0', '0.2', '0.3', '0.4', '0.6', '0.8', '1.0'}, ...
    {'0', '0.25', '0.5', '0.75', '1.0', '1.25', '1.5'}, ...
    {'0', '0.25', '0.5', '0.75', '1.0', '1.25', '1.5'}, ...
    {'-0.5', '-0.25', '0', '0.25', '0.5', '0.75', '1.0'}, ...
    [], ...
    [], ...
    [], ...
    [], ...
    [], ...
    {'50', '60', '70', '80', '90', '100'}, ...
    {'50', '60', '70', '80', '90', '100'}, ...
    [], ...
    {'0', '0.2', '0.3', '0.4', '0.6', '0.8', '1.0'}, ...
    };
%
%
n_cols_figure   = 2;
journal         = 'mrm';
res_dpi         = 300;
%
%
m_upper = 0.2;
ph      = 1;         % plot
phm     = 0.5;
m_lower = 0.5;
%
m_left  = 0.5;
pw      = 1;
pwm     = 0.5;
m_right = 0.2;
%
axis_factor     = 1;

% Axis setup
fh = m_upper    + n_row * ph + (n_row - 1) * phm + m_lower;
%
fw = m_left     + n_col * pw + (n_col - 1) * pwm + m_right;
%
m_upper = m_upper / fh;
ph      = ph / fh;
phm     = phm / fh;
m_lower = m_lower / fh;
%
m_left  = m_left / fw;
pw      = pw / fw;
pwm    = pwm / fw;
m_right = m_right / fw;
%
imsize     = my_imsize(n_cols_figure, fh / fw, journal);
%
fs_axis         = 6.5;
fs_lbl          = 6.5;
%
lw_axis         = 0.8;
lw_gt           = 0.8;
col_gt          = [0.7 0.4 0.4];
tl              = [0.03 0.03];
%
p2_gt           = [0.87 0.14 0 0]; % OD = 0.05, 0.5, 1, 1
%
m_bs.ms     = 3.5;
m_bs.mec    = 'k';
x           = [1 2 3 4];


%%%------------------------------------------------------------------
% SIMULATION
%%%------------------------------------------------------------------


output_fn = fullfile(output_dir, 'smr_fig4_boot.mat');
if (~exist(output_fn, 'file') || do_overwrite)
    boot_mp       = zeros(n_mp, opt_boot.n_boot, 14); % msr, then add p2
    boot_s        = zeros(n_mp, opt_boot.n_boot, xps.n);
    boot_s_fit    = boot_s;
    %
    for c_mp = 1:n_mp
        
        
        disp(['mp = ' mp_list{c_mp}]);
        disp(' ');
        %
        mpw = smr_mp_watson(mp_list{c_mp});
        %
        for c = 1:2
            if (c == 1)
                opt_boot.opt_fit.init_guess = [];
                [boot_mp1, ~, ~]                 = smr_fig4_boot(xps, mpw, opt_boot);
            else
                opt_boot.opt_fit.init_guess = median(boot_mp1(:,1:12)) .* [1 1 1e9 1e9 1 1 1 1 1 1 1e3 1e3];
                [boot_mp2, boot_s2, boot_s_fit2] = smr_fig4_boot(xps, mpw, opt_boot);
            end
        end
        %
        p2_boot = smr_p2m2p2(boot_mp2(:,6:10));
        %
        boot_mp   (c_mp, :, :) = [boot_mp2(:,1:13) p2_boot];
        boot_s    (c_mp, :, :) = boot_s2;
        boot_s_fit(c_mp, :, :) = boot_s_fit2;
    end
    
    % Save
    %
    save(output_fn, ...
        'boot_mp')
    disp(['Wrote: ' output_fn]);    
end
my_struct = load(output_fn);

%%%------------------------------------------------------------------
% PLOTTING
%%%------------------------------------------------------------------


%%% Figure
%
figure(675)
clf
set(gcf, 'color', 'w')
%
%
pars_sd     = zeros(n_par, n_mp);
c_p         = 1;
%
for c_row = 1:n_row
    for c_col = 1:n_col
        
        if (c_p > n_par), continue; end
        
        c_par   = c_par_list(c_p);
        c_par_w = c_par_list_w(c_p);
        
        %%% Bar plot
        %
        ax_l = m_left + (c_col - 1) * (pw + pwm);
        ax_b = 1 - m_upper - ph - (c_row - 1) * (ph + phm);
        ax_w = pw * axis_factor;
        ax_h = ph * axis_factor;
        af_adjust = - pw * (axis_factor - 1) / 2; % Correct for axis factor
        ax_l = ax_l + af_adjust;
        ax_b = ax_b + af_adjust;
        axes('position', [ax_l ax_b ax_w ax_h]);
        
        
        % Obtain gt, means and std
        vals    = zeros(n_mp, 200);
        vals_gt = zeros(n_mp, 1);
        vals_m  = vals_gt;
        vals_sd = vals_gt;
        %
        %
        for c_mp = 1:n_mp
            %
            mpw     = smr_mp_watson(mp_list{c_mp});
            if (c_p == n_par)
                vals_gt(c_mp) = p2_gt(c_mp);
            else
                vals_gt(c_mp) =  mpw(c_par_w) * par_sc(c_par);
            end            
            %
            tmp = squeeze(my_struct.boot_mp(c_mp, :, c_par)) * par_sc(c_par);
            %
            vals   (c_mp, :) = tmp;
            vals_m (c_mp) = mean(tmp);
            vals_sd(c_mp) = std(tmp);
            %
            pars_sd(c_p, c_mp) = vals_sd(c_mp);
            %
        end
        
        % Plot ground truth
        for c_mp = 1:n_mp
            plot(linspace(c_mp-0.5, c_mp+0.5, 10), repmat(vals_gt(c_mp), [10 1]), '--', ...
                'color', col_gt, 'linewidth', lw_gt)
            hold on
        end        
        
        % Plot beeswarms
        for c_mp = 1:n_mp
            my_plot_beeswarm(x(c_mp), {vals(c_mp,:)}, [], [], m_bs, 0);
        end
        
        % Axis and labels
        set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', fs_axis, 'linewidth', lw_axis, ...
            'fontweight', 'bold', 'ticklength', tl)
        axis square
        ylabel(par_list{c_par}, 'fontsize', fs_lbl)
        xlim([0.5 n_mp+0.5])
        ylim(par_lims{c_par})
        yl = ylim;
        %
        set(gca, 'xtick', 1:n_mp, 'xticklabel', xticklbl, 'xticklabelrot', 45, ...
            'fontsize', fs_axis, 'fontweight', 'bold');
        set(gca, 'ytick', linspace(yl(1), yl(2), numel(par_tick_lbl{c_par})), 'yticklabel', par_tick_lbl{c_par}, ...
            'fontsize', fs_axis, 'fontweight', 'bold');
        
        %
        c_p = c_p + 1;
        
    end
end


%%% PRINT FIGURE TO FILE
output_filename = fullfile(output_dir, output_name);
save_current_fig_to_file(output_name, output_dir, imsize, res_dpi);
disp(['Wrote: ' output_filename]);
system(['open ' output_filename]);

end