function smr_fig3(output_dir)

if (nargin < 1), output_dir = fullfile(pwd, 'figures', 'output'); end

do_overwrite            = 0;
output_name             = 'smr_fig3.tiff';

%%%------------------------------------------------------------------
% PREPARE
%%%------------------------------------------------------------------


%%% EP and MP
%
ep_list = {...
    'opt_10_shell', ...
    'in_vivo', ...
    'opt_lte', ...
    'lampinen_hbm_2019'};
n_ep = numel(ep_list);
%
mp      = smr_optimize_mp('b');

%%% NRV iteration
%
opt_optimize            = smr_optimize_opt_derive(smr_optimize_opt);
snr                     = opt_optimize.snr;
%
opt_fit                 = smr_opt;
opt_fit.n_rep           = 2;
%
opt_iter.n_iter         = 10;
opt_iter.snr            = snr;
opt_iter.opt_fit        = opt_fit;
step = 0.0125;
opt_iter.fs             = (step/2):step:(1-step/2);
%

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
b_xlbl = 'fixed \itf\rm\bf_S (in fitting)';
b_ylbl = 'NRV';
%
n_cols_figure   = 1;
journal         = 'mrm';
res_dpi         = 300;
%
%
m_left      = 0.3;
pw          = 1;
m_right     = 0.8;
%
m_upper     = 0.15;
ph          = 1;
m_lower     = 0.3;
%
axis_factor     = 1;

% Axis setup
fh = m_upper    + ph + m_lower;
fw = m_left     + pw + m_right;
%
%
m_upper = m_upper / fh;
ph      = ph / fh;
m_lower = m_lower / fh;
%
m_left  = m_left / fw;
pw      = pw / fw;
m_right = m_right / fw;
%
imsize     = my_imsize(n_cols_figure, fh / fw, journal);
%
fs_legend       = 7;
fs_axis         = 8;
fs_lbl          = 9;
%
nrv_min = 0.9;
nrv_max = 2;
lw_gt = 0.8;
col_gt = 0.5 + [0 0 0];
%
lw_axis         = 0.8;
lw_plot_list    = 1.2 * [1 1 1 1];
%
yticklbl        = {'1.0',  '1.2' , '1.4' , '1.6' , '1.8' , '2.0'};
xticklbl        = {'0.00', '0.25', '0.50', '0.75', '1.00'};


%%%------------------------------------------------------------------
% NRV-ITERATION
%%%------------------------------------------------------------------

for c_ep = 1:n_ep
    
    output_fn = fullfile(output_dir, ['smr_fig3_' ep_list{c_ep} '.mat']);
    if (exist(output_fn, 'file') && ~do_overwrite)
        disp(['Exists: ' output_fn]);
        continue;
    end
    %
    disp(['ep = ' ep_list{c_ep}]);
    disp(' ');
    
    % Adjust SNR
    xps = my_ep2xps(smr_ep(ep_list{c_ep}));
    tacq_factor = smr_optimize_xps2tacq(xps, opt_optimize) / opt_optimize.tacq_limit;
    opt_iter.snr = opt_optimize.snr / sqrt(tacq_factor);
    
    % Iter
    %
    [ssr_vals, mp_vals] = smr_fig3_iterate(xps, mp, opt_iter);
    
    % Save
    %
    save(output_fn, 'ssr_vals', 'mp_vals', 'xps')
    disp(['Wrote: ' output_fn]);
end


%%%------------------------------------------------------------------
% PLOTTING
%%%------------------------------------------------------------------


%%% Figure
%
figure(897)
clf
set(gcf, 'color', 'w')

%
%
ax_l = m_left;
ax_b = 1 - m_upper - ph;
ax_w = pw * axis_factor;
ax_h = ph * axis_factor;
af_adjust = - pw * (axis_factor - 1) / 2; % Correct for axis factor
ax_l = ax_l + af_adjust;
ax_b = ax_b + af_adjust;
axes('position', [ax_l ax_b ax_w ax_h]);
%
%


%%% Legend
for c_ep = 1:n_ep
    plot(-1, -1, '-', ...
        'color', ep_col_list{c_ep_legend_list(c_ep)}, ...
        'linewidth', lw_plot_list(c_ep_legend_list(c_ep)));
    hold on
end


%%% Plot straight line (GT)
plot(repmat(mp(2), [1 10]), linspace(nrv_min, nrv_max, 10), '--', 'color', col_gt, 'linewidth', lw_gt)
hold on


%%% Plot NRV
c_ep_plot_list = [4 3 1 2];
for c_ep = 1:n_ep
    
    c_ep_plot = c_ep_plot_list(c_ep);
    
    % XPS
    xps = my_ep2xps(smr_ep(ep_list{c_ep_plot}));
    tacq_factor = smr_optimize_xps2tacq(xps, opt_optimize) / opt_optimize.tacq_limit;
    opt_iter.snr = opt_optimize.snr / sqrt(tacq_factor);
    
    % Load
    output_fn = fullfile(output_dir, ['smr_fig3_' ep_list{c_ep_plot} '.mat']);
    my_struct = load(output_fn);
    ssr_vals  = my_struct.ssr_vals;
    
    % Calculate NRV
    sigma_noise = 1 / opt_iter.snr;
    nrv         = (ssr_vals / (xps.n - 11)) / sigma_noise^2;
    nrv         = nrv ./ repmat(min(nrv, [], 2), [1 80]);
    
    % Obtain a smoothed average line
    n_sd = 1;
    c_ol = abs(nrv - repmat(median(nrv), [opt_iter.n_iter 1])) > repmat(n_sd * std(nrv), [opt_iter.n_iter 1]);
    tmp = nrv;
    tmp(c_ol) = NaN;
    nrv_m   = smooth(smooth(nanmedian(tmp)));
    
    % Plot
    plot(opt_iter.fs, nrv_m, '-', ...
        'color', ep_col_list{c_ep_plot}, ...
        'linewidth', lw_plot_list(c_ep_plot));
    hold on
end

%
l = legend(name_list(c_ep_legend_list), 'box', 'off', 'fontsize', fs_legend, 'fontweight', 'normal');
pos = get(l, 'pos');
set(l, 'pos', [1.6 * pos(1) 0.35 * pos(2) pos(3:4)]);



set(gca, 'xtick', linspace(0,1,5), 'xticklabel', {'0', '0.25', '0.5', '0.75', '1'})
set(gca, 'box', 'off', 'tickdir', 'out');
set(gca, 'xticklabel', xticklbl, 'yticklabel', yticklbl)
axis square
ylabel(b_ylbl, 'fontsize', fs_lbl, 'fontweight', 'bold')
xlabel(b_xlbl, 'fontsize', fs_lbl, 'fontweight', 'bold')
set(gca, 'box', 'off', 'tickdir', 'out', 'fontweight', 'bold', 'linewidth', lw_axis, 'fontsize', fs_axis)
axis([0 1 nrv_min nrv_max])

%
%%% PRINT FIGURE TO FILE
output_filename = fullfile(output_dir, output_name);
save_current_fig_to_file(output_name, output_dir, imsize, res_dpi);
disp(['Wrote: ' output_filename]);
system(['open ' output_filename]);

end
