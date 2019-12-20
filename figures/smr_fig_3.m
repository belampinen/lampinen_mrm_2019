function smr_fig_3(output_dir)

if (nargin < 1), output_dir = fullfile(pwd, 'figures', 'output'); end

output_name             = 'smr_fig_3.tiff';


%%%------------------------------------------------------------------
% PREPARE
%%%------------------------------------------------------------------


%%% MISC
ep_list = {...
    'opt_10_shell', ...
    'in_vivo', ...
    'opt_lte', ...
    'lampinen_hbm_2019'};
do = [1 0 1 0];
n_ep = numel(ep_list);
%
mp_list = {'a', 'b', 'c'};
n_mp    = numel(mp_list);
%
opt_optimize        = smr_optimize_opt_derive(smr_optimize_opt);
%
n_iter              = 10;
step                = 0.0125;
fs_fix              = (step/2):step:(1-step/2);


%%% Plotting
%
ep_col_list   = {[0.7 0.4 0.4], [0 0 0], [0.4 0.7 0.4], [0.4 0.4 0.7]};
c_ep_legend_list = [1 2 3 4];
%
name_list   = {...
    'CRLB-optimized', ...
    'In vivo protocol', ...
    'CRLB-optimized (LTE only)', ...
    'Multiple TE only for low \itb'};
%
b_xlbl    = 'fixed \itf\rm\bf_S (in fitting)';
title_list = {'(A) Prior set A', '(B) Prior set B', '(C) Prior set C'};
y_lbl     = 'NRV';
%
txt_start = '\bftrue \itf\rm\bf_S';
txt_str = {...
    [txt_start ' (0.45)'], ...
    [txt_start ' (0.15)'], ...
    [txt_start ' (0.40)'], ...
    };
%
n_cols_figure   = 2;
journal         = 'mrm';
res_dpi         = 600;
%
%
m_left      = 0.3;
th          = 0.2;
thm         = 0.1;
pw          = 1;
pwm         = 0.4;
m_right     = 0.8;
%
m_upper     = 0.06;
ph          = 1;
% phm         = 0.5;
m_lower     = 0.26;
%
tm          = 0.25;
tw          = 2; % A, B, C, D...
%
axis_factor     = 1;
%
leg_start_l     = 0.79;
leg_start_b     = 0.35;

% Axis setup
fh = m_upper    + th + thm + ph + m_lower;
fw = m_left     + n_mp * pw + (n_mp - 1) * pwm + m_right;
%
%
m_upper = m_upper / fh;
th      = th / fh;
thm     = thm / fh;
ph      = ph / fh;
m_lower = m_lower / fh;
%
m_left  = m_left / fw;
pw      = pw / fw;
pwm     = pwm / fw;
m_right = m_right / fw;
%
tm      = tm / fw;
tw      = tw / fw;
%
imsize     = my_imsize(n_cols_figure, fh / fw, journal);
%
fs_ref          = 5;
fs_legend       = fs_ref + 1;
fs_axis         = fs_ref + 1;
fs_lbl          = fs_ref + 2;
fs_title        = fs_ref + 3;
fs_text         = fs_ref + 1;
%
nrv_min = 0.9;
nrv_max = 2;
% nrv_max = 4;
lw_gt   = 0.8;
col_gt  = 0.4 + [0 0 0];
%
lw_axis         = 0.8;
lw_plot_list    = 1.2 * [1 1 1 1];
%
yticklbl        = {'1.0',  '1.2' , '1.4' , '1.6' , '1.8' , '2.0'};
xticklbl        = {'0.00', '0.25', '0.50', '0.75', '1.00'};
%    
txt_offset_b = 0.07;
txt_offset_l = -0.16;



%%%------------------------------------------------------------------
% PERFORM
%%%------------------------------------------------------------------


%%% Figure
%
figure(897)
clf
set(gcf, 'color', 'w')
%


for c_mp = 1:n_mp
    
    mp      = smr_mp(mp_list{c_mp});
    %
    l_pos   = m_left + (c_mp - 1) * (pw + pwm);
    
    %%% Title
    %
    t_l = max(l_pos - tm, 0);
    t_b = 1 - m_upper - th;
    t_w = tw;
    t_h = th;
    annotation(...
        'textbox', [t_l t_b t_w t_h], ...
        'string', title_list{c_mp}, ...
        'fontweight', 'normal', ...
        'FitBoxToText', 'on', ...
        'linestyle', 'none', ...
        'horizontalalignment', 'left', ...
        'verticalalignment', 'middle', ...
        'fontsize', fs_title);    
    
    
    %%% Axes
    %
    ax_l = l_pos;
    ax_b = 1 - m_upper - th - thm - ph;
    ax_w = pw * axis_factor;
    ax_h = ph * axis_factor;
    af_adjust = - pw * (axis_factor - 1) / 2; % Correct for axis factor
    ax_l = ax_l + af_adjust;
    ax_b = ax_b + af_adjust;
    axes('position', [ax_l ax_b ax_w ax_h]);
    %
    
    if (c_mp == n_mp)
        %%% Prepare legend
        for c_ep = 1:n_ep
            plot(-1, -1, '-', ...
                'color', ep_col_list{c_ep_legend_list(c_ep)}, ...
                'linewidth', lw_plot_list(c_ep_legend_list(c_ep)));
            hold on
        end
    end
  
    %%% Plot straight line (GT)
    plot(repmat(mp(2), [1 10]), linspace(nrv_min, nrv_max, 10), '--', 'color', col_gt, 'linewidth', lw_gt)
    hold on
    text(mp(2) + txt_offset_l, nrv_max + txt_offset_b, ...
       txt_str{c_mp}, ...
        'fontsize', fs_text, 'color', col_gt);
    
    
    %%% Plot NRV
    c_ep_plot_list = [4 3 1 2];
    for c_ep = 1:n_ep
        
        c_ep_plot = c_ep_plot_list(c_ep);
        
        % XPS
         %
        if (do(c_ep))
            ep = smr_ep(ep_list{c_ep});
        else
            ep = smr_ep(ep_list{c_ep});
        end  
        xps         = my_ep2xps(ep); 
        tacq_factor = smr_optimize_xps2tacq(xps, opt_optimize) / opt_optimize.tacq_limit;
        snr_ep      = opt_optimize.snr / sqrt(tacq_factor);
        
        % Load
        output_fn = fullfile(output_dir, ['smr_fig_3_sim_' mp_list{c_mp} '_' ep_list{c_ep_plot} '.mat']);
        my_struct = load(output_fn);
        ssr_vals  = my_struct.ssr_vals;
        
        % Calculate NRV
        sigma_noise = 1 / snr_ep;
        nrv         = (ssr_vals / (xps.n - 11)) / sigma_noise^2;
        nrv         = nrv ./ repmat(min(nrv, [], 2), [1 80]);
        
        % Obtain a smoothed average line
        n_sd = 1;
        c_ol = abs(nrv - repmat(median(nrv), [n_iter 1])) > repmat(n_sd * std(nrv), [n_iter 1]);
        tmp = nrv;
        tmp(c_ol) = NaN;
        nrv_m   = smooth(smooth(nanmedian(tmp)));
        
        % Plot
        plot(fs_fix, nrv_m, '-', ...
            'color', ep_col_list{c_ep_plot}, ...
            'linewidth', lw_plot_list(c_ep_plot));
        hold on
    end
    
    %
    if (c_mp == n_mp)
        l = legend(name_list(c_ep_legend_list), 'box', 'off', ...
            'fontsize', fs_legend, 'fontweight', 'normal', 'location', 'southeast');
        pos = get(l, 'pos');
        set(l, 'pos', [leg_start_l leg_start_b pos(3:4)]);
        set(l, 'ItemTokenSize', [10 18]);
    end
    
    set(gca, 'xtick', linspace(0,1,5), 'xticklabel', {'0', '0.25', '0.5', '0.75', '1'})
    set(gca, 'box', 'off', 'tickdir', 'out');
    set(gca, 'xticklabel', xticklbl, 'yticklabel', yticklbl)
    axis square
    ylabel(y_lbl, 'fontsize', fs_lbl, 'fontweight', 'bold')
    xlabel(b_xlbl, 'fontsize', fs_lbl, 'fontweight', 'bold')
    set(gca, 'box', 'off', 'tickdir', 'out', 'fontweight', 'bold', 'linewidth', lw_axis, 'fontsize', fs_axis)
    axis([0 1 nrv_min nrv_max])
end

%
%%% PRINT FIGURE TO FILE
output_filename = fullfile(output_dir, output_name);
save_current_fig_to_file(output_name, output_dir, imsize, res_dpi);
disp(['Wrote: ' output_filename]);
system(['open ' output_filename]);

end
