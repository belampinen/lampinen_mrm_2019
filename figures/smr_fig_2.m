function smr_fig_2(output_dir)

if (nargin < 1), output_dir = fullfile(pwd, 'figures', 'output'); end

output_name = 'smr_fig_2.tiff';

%%%------------------------------------------------------------------
% PREPARE
%%%------------------------------------------------------------------

opt = smr_optimize_opt_derive(smr_optimize_opt);
%
ep_list_shell = {...
    'opt_4_shell', ...
    'opt_5_shell', ...
    'opt_6_shell', ...
    'opt_7_shell', ...
    'opt_8_shell', ...
    'opt_9_shell', ...
    'opt_10_shell'};
%
ep_list_shell_name = {...
    '4', ...
    '5', ...
    '6', ...
    '7', ...
    '8', ...
    '9', ...
    '10'};
%
ep_list_bdelta = {...
    'opt_lte', ...
    'opt_dde', ...
    'opt_ste', ...
    'opt_10_shell', ...
    };
%
ep_list_bdelta_name = {...
    '1', ...
    '1/-0.5', ...
    '1/0', ...
    'any'};

b_delta_lbl = '\bf\itb\rm\bf_{\Delta} combinations';
rot_bdelta  = 0;
%
%
gmax_list  = [75 200 Inf];
gmax_adder = {'_75', '_200', '_inf'};
gmax_leg_str     = {...
    '\itg\rm\bf_{max} = 80 mT/m', ...
    '\itg\rm\bf_{max} = 200 mT/m', ...
    '\itg\rm\bf_{max} = Inf', ...
    };
n_gmax     = numel(gmax_list);
%
ep_list_bmax = {...
    'opt_b2000', ...
    'opt_b2500', ...
    'opt_b3000', ...
    'opt_b4000', ...
    'opt_b5000', ...
    'opt_b6000', ...
    'opt_b7000'};
bmax_list = [2 2.5 3 4 5 6 7];
bmax_tick = [2 3 4 5 6 7];
%
%
ep_list_temin = {...
    'opt_te10', ...
    'opt_te20', ...
    'opt_te30', ...
    'opt_te40', ...
    'opt_te50', ...
    'opt_te60', ...
    'opt_te70', ...
    'opt_te80', ...
    'opt_te90', ...
    'opt_te100', ...
    };
temin_list = [10 20 30 40 50 60 70 80 90 100];
temin_tick = [10 25 40 55 70 85 100];

%
%
%
a_str           = '(A) - Maximal b-value';
b_str           = '(B) - Shapes of the b-tensor';
c_str           = '(C) - Minimal echo time';
d_str           = '(D) - Number of shells';%
sdw_str = 'V_{W}';
%
%
shell_lbl   = 'Number of shells';
bmax_lbl    = ['\itb\rm\bf_{max} [ms/' char(181) 'm^2]'];
temin_lbl   = 'TE_{min} [ms]';
%
ymax_vw     = 2.2;
ar_hw       = ymax_vw / 1.7;
ytick       = [0 0.5 1 1.5];
yticklbl    = {'0.0', '0.5', '1.0', '1.5'};
rep_4sh     = 2.0;
rep_5sh     = 1.75;
rep_lte     = 1.85;
%
%
yl_pos = 1.1; 
%
l_sc = 0.14;
b_sc = 0.85;
%
l_sc_b = 0.915;
b_sc_b = 2.65;
w_sc_b = 0.3;
h_sc_b = 0.3;


%%%------------------------------------------------------------------
% PERFORM
%%%------------------------------------------------------------------

n_cols_figure   = 2;
journal         = 'mrm';
res_dpi         = 600;
%
%
%
m_left   = 0.5;
pw       = 1;
pwm      = 0.6;
m_right  = 0.2;
%
m_upper = 0.1;
th      = 0.2;
thm     = 0; % 0.2;  % for broken bars
ph      = 1 * ar_hw;
m_lower = 0.45;
%
tm      = 0.25;
tw      = 2; % A, B, C, D...
%
axis_factor     = 1;
%

% Axis setup
fh = m_upper    + th + thm + ph + m_lower;
fw = m_left     + 4 * pw + 3 * pwm + m_right;
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
tm      = tm / fw;
tw      = tw / fw;
%
imsize     = my_imsize(n_cols_figure, fh / fw, journal);
%
fs_title        = 9;
fs_other        = 7;
fs_legend_bm    = fs_other - 1;
fs_axis         = fs_other;
fs_lbl          = fs_other;
fs_cb           = fs_other - 1;
%
lw_axis         = 1;
lw_bar          = 0.5;
%
ntick_cb = 3;
lw_cb = 0.5;
tl_cb = 0.025;
tl_axis = 0.025 + [0 0];
%
bar_color       = 0.8 + [0 0 0]; %[0.2 0.8 1];
bw1             = 0.3;
bw2             = 0.2;

%
% CRLB-plots
lw_crlb          = lw_bar;
%
gmax_marker_list = {'o', 's', 'd'};
gmax_ms_list     = 4 + [0 0 0];
gmax_ms_list_col = 3 + [0 0 0];
lw_bm_marker     = 0.2;


rot_bmax   = 0;
rot_temin  = 0;


%%% Figure
%
figure(967)
clf
set(gcf, 'color', 'w')


%%% A: B-value
%
%
l_pos = m_left;
%
% A annotation
t_l = l_pos - tm;
t_b = 1 - m_upper - th;
t_w = tw;
t_h = th;
annotation(...
    'textbox', [t_l t_b t_w t_h], ...
    'string', a_str, ...
    'fontweight', 'normal', ...
    'FitBoxToText', 'on', ...
    'linestyle', 'none', ...
    'horizontalalignment', 'left', ...
    'verticalalignment', 'middle', ...
    'fontsize', fs_title);
%
%
% A axes
ax_l = l_pos;
ax_b = 1 - m_upper - th - thm - ph;
ax_w = pw * axis_factor;
ax_h = ph * axis_factor;
af_adjust = - pw * (axis_factor - 1) / 2; % Correct for axis factor
ax_l = ax_l + af_adjust;
ax_b = ax_b + af_adjust;
axes('position', [ax_l ax_b ax_w ax_h]);
l_pos = l_pos + pw + pwm;

% B-values
% Obtain values
n_ep = numel(ep_list_bmax);
%
vw = zeros(n_ep, n_gmax);
for c_gmax = 1:n_gmax
    for c_ep = 1:n_ep
        ep                  = smr_ep([ep_list_bmax{c_ep} gmax_adder{c_gmax}]);
        vw(c_ep,c_gmax)     = smr_optimize_ep2vw(ep, opt);
    end
end

vw_ch = 100 * (vw(end,2) - vw(end,1)) / vw(end,1);
disp(['80 --> 200 mT/m reuduced Vw ~ ' num2str(vw_ch) ' percent']);

% Create legend (sorted by shell #)
for c_gmax = 1:n_gmax
    plot(zeros(100,1)-1, ones(100,1), [gmax_marker_list{c_gmax} 'k'], ...
        'markersize', gmax_ms_list(c_gmax), 'linewidth', lw_bm_marker, ...
        'markerfacecolor', 'k');
    hold on
end

l = legend(gmax_leg_str, 'box', 'off', ...
    'fontsize', fs_legend_bm);
pos = get(l, 'pos');
set(l, 'pos', [(ax_l + l_sc * pw) b_sc * pos(2) pos(3:4)]) % pos(3) pos(4)]) % lbwh


% Plot
for c_gmax = 1:n_gmax
    plot(bmax_list, vw(:,c_gmax), '-k', 'linewidth', lw_crlb)
    hold on
    plot(bmax_list, vw(:,c_gmax), [gmax_marker_list{c_gmax} 'k'], 'markersize', gmax_ms_list(c_gmax), ...
        'linewidth', lw_bm_marker, 'markerfacecolor', 'k');
end


% Axis and labels
set(gca, 'xtick', bmax_tick, 'xticklabel', bmax_tick, ...
    'box', 'off', 'tickdir', 'out', 'fontweight', 'bold', 'linewidth',...
    lw_axis, 'fontsize', fs_axis, 'ticklength', tl_axis, ...
    'ytick', ytick, 'yticklabel', yticklbl, ...
    'xticklabelrot', rot_bmax)
%     end
xlim([bmax_list(1) bmax_list(end)])
ylim([0 ymax_vw])
ylabel(sdw_str, 'fontsize', fs_lbl)
xlabel(bmax_lbl, 'fontsize', fs_lbl)
%
yl = get(gca, 'ylabel');
p  = get(yl, 'position');
set(yl, 'position', [p(1) yl_pos p(3)])


%%% B: Shapes of the b-tensor
%
% B annotation
t_l = l_pos - tm;
t_b = 1 - m_upper - th;
t_w = tw;
t_h = th;
annotation(...
    'textbox', [t_l t_b t_w t_h], ...
    'string', b_str, ...
    'fontweight', 'normal', ...
    'FitBoxToText', 'on', ...
    'linestyle', 'none', ...
    'horizontalalignment', 'left', ...
    'verticalalignment', 'middle', ...
    'fontsize', fs_title);
%
% B axes
ax_l = l_pos;
ax_b = 1 - m_upper - th - thm - ph;
ax_w = pw * axis_factor;
ax_h = ph * axis_factor;
af_adjust = - pw * (axis_factor - 1) / 2; % Correct for axis factor
ax_l = ax_l + af_adjust;
ax_b = ax_b + af_adjust;
axes('position', [ax_l ax_b ax_w ax_h]);
l_pos = l_pos + pw + pwm;

% B-tensor shapes
%
% Obtain values
n_ep = numel(ep_list_bdelta);
vw   = zeros(n_ep, 1);
for c_ep = 1:n_ep
    ep          = smr_ep(ep_list_bdelta{c_ep});
    vw(c_ep)    = smr_optimize_ep2vw(ep, opt);
end
disp(['LTE-only: Vw = ' num2str(vw(1))]);
%
%
vw_ch = 100 * (vw(2) - vw(1)) / vw(1);
disp(['lte --> dde shells reuduced Vw ~ ' num2str(vw_ch) ' percent']);
%
vw(1) = rep_lte; % Replace display value to fit in plot


% Plot
bar(1:n_ep, vw, 'facecolor', bar_color, 'linewidth', lw_bar, 'barwidth', bw2)

% Axis and labels
set(gca, 'xtick', 1:n_ep, 'xticklabel', ep_list_bdelta_name, ...
    'box', 'off', 'tickdir', 'out', 'fontweight', 'bold', 'linewidth', lw_axis, 'fontsize', fs_axis, ...
    'xticklabelrot', rot_bdelta, 'ticklength', tl_axis, ...
     'ytick', ytick, 'yticklabel', yticklbl)
xlim([0 n_ep+1])
ylim([0 ymax_vw])
if (~isempty(b_delta_lbl))
    xlabel(b_delta_lbl, 'fontsize', fs_lbl)
end
ylabel(sdw_str, 'fontsize', fs_lbl)
%
yl = get(gca, 'ylabel');
p  = get(yl, 'position');
set(yl, 'position', [p(1) yl_pos p(3)])




%%% C: TEmin
%
% C annotation
t_l = l_pos - tm;
t_b = 1 - m_upper - th;
t_w = tw;
t_h = th;
annotation(...
    'textbox', [t_l t_b t_w t_h], ...
    'string', c_str, ...
    'fontweight', 'normal', ...
    'FitBoxToText', 'on', ...
    'linestyle', 'none', ...
    'horizontalalignment', 'left', ...
    'verticalalignment', 'middle', ...
    'fontsize', fs_title);
%
% C axes
ax_l = l_pos;
ax_b = 1 - m_upper - th - thm - ph;
ax_w = pw * axis_factor;
ax_h = ph * axis_factor;
af_adjust = - pw * (axis_factor - 1) / 2; % Correct for axis factor
ax_l = ax_l + af_adjust;
ax_b = ax_b + af_adjust;
axes('position', [ax_l ax_b ax_w ax_h]);
l_pos = l_pos + pw + pwm;


%%% TEmin
% Obtain values
n_ep = numel(ep_list_temin);
%
vw = zeros(n_ep, n_gmax);
for c_gmax = 1:n_gmax
    for c_ep = 1:n_ep
        ep                  = smr_ep([ep_list_temin{c_ep} gmax_adder{c_gmax}]);
        vw(c_ep,c_gmax)     = smr_optimize_ep2vw(ep, opt);
    end
end
%
%
%
temin_used  = zeros(n_gmax, n_ep);
b_max_used  = temin_used;
b_max_temin = zeros(n_gmax, n_ep);
%
c_gmax_list = [3 2 1];
for c_g = 1:3
    for c_ep = 1:n_ep
        %
        ep     = smr_ep([ep_list_temin{c_ep} gmax_adder{c_gmax_list(c_g)}]);
        
        %
        temin_used (c_gmax_list(c_g), c_ep) = min(ep(:,3));
        b_max_used (c_gmax_list(c_g), c_ep) = max(ep(:,1));
        %
        b_max_temin(c_gmax_list(c_g), c_ep) = max(ep( ep(:,3) == temin_used(c_gmax_list(c_g), c_ep), 1));
    end
end
%
n_b = 200;
b_vec = linspace(0, 7, n_b);
cm    = ff_cmap_viridis(n_b + 30);
cm    = cm(11:end,:);


% Create legend (sorted by shell #)
for c_gmax = 1:n_gmax
    p = plot(zeros(100,1)-1, ones(100,1), [gmax_marker_list{c_gmax} 'k'], ...
        'markersize', gmax_ms_list(c_gmax), 'linewidth', lw_bm_marker, ...
        'markerfacecolor', 'k');
    hold on
end
l = legend(gmax_leg_str, 'box', 'off', ...
    'fontsize', fs_legend_bm);
pos = get(l, 'pos');

set(l, 'pos', [(ax_l + l_sc * pw) b_sc * pos(2) pos(3:4)]) % pos(3) pos(4)]) % lbwh


%
c_gmax_list = [3 2 1];
for c_g  = 1:3
    
    c_gmax = c_gmax_list(c_g);
    
    % Plot
    plot(temin_list, vw(:,c_gmax), 'k-', 'linewidth', lw_crlb)
    hold on        
    
    % Points with color on black
    for c_ep = 1:n_ep
        % black
        plot(temin_list(c_ep), vw(c_ep,c_gmax), [gmax_marker_list{c_gmax} 'k'], ...
            'markersize', gmax_ms_list(c_gmax));
        % color
        col_ind = find(abs(b_max_temin(c_gmax, c_ep) - b_vec) == min(abs(b_max_temin(c_gmax, c_ep) - b_vec)));
        if (numel(col_ind) > 1), col_ind = col_ind(1); end
        %
        plot(temin_list(c_ep), vw(c_ep,c_gmax), gmax_marker_list{c_gmax}, ...
            'markerfacecolor', cm(col_ind, :), 'markeredgecolor',  cm(col_ind, :), ...
            'markersize', gmax_ms_list_col(c_gmax));
    end
end


% Axis and labels
set(gca, 'xtick', temin_tick, 'xticklabel', temin_tick, ...
    'box', 'off', 'tickdir', 'out', 'fontweight', 'bold', ...
    'linewidth', lw_axis, 'fontsize', fs_axis, 'xticklabelrot', rot_temin, ...
    'ticklength', tl_axis,  'ytick', ytick, 'yticklabel', yticklbl)
%     end
xlim([temin_list(1) temin_list(end)])
ylim([0 ymax_vw])
ylabel(sdw_str, 'fontsize', fs_lbl)
xlabel(temin_lbl, 'fontsize', fs_lbl)
%
yl = get(gca, 'ylabel');
p  = get(yl, 'position');
set(yl, 'position', [p(1) yl_pos p(3)])

set(gcf, 'colormap', cm)
cb = colorbar;
set(cb, 'fontsize', fs_cb, 'fontweight', 'bold', 'tickdir', 'out', ...
    'ticks', linspace(0,1,ntick_cb), 'ticklabels', {'0.0', '3.5', '7.0'}, ...
    'linewidth', lw_cb, 'ticklength', tl_cb)
%
p1 = get(cb, 'pos');
%
set(cb, 'pos', [l_sc_b * p1(1) b_sc_b * p1(2) w_sc_b * p1(3) h_sc_b * p1(4)])


%%% Panel D: Combinations
%
% D annotation
t_l = l_pos - tm;
t_b = 1 - m_upper - th;
t_w = tw;
t_h = th;
annotation(...
    'textbox', [t_l t_b t_w t_h], ...
    'string', d_str, ...
    'fontweight', 'normal', ...
    'FitBoxToText', 'on', ...
    'linestyle', 'none', ...
    'horizontalalignment', 'left', ...
    'verticalalignment', 'middle', ...
    'fontsize', fs_title);
%
% D axes
ax_l = l_pos;
ax_b = 1 - m_upper - th - thm - ph;
ax_w = pw * axis_factor;
ax_h = ph * axis_factor;
af_adjust = - pw * (axis_factor - 1) / 2; % Correct for axis factor
ax_l = ax_l + af_adjust;
ax_b = ax_b + af_adjust;
axes('position', [ax_l ax_b ax_w ax_h]);


%%% Combinations
%
% Obtain values
n_ep = numel(ep_list_shell);
vw   = zeros(n_ep, 1);
for c_ep = 1:n_ep
    ep          = smr_ep(ep_list_shell{c_ep});
    vw(c_ep)    = smr_optimize_ep2vw(ep, opt);
end

disp(['4 shells: Vw = ' num2str(vw(1))]);
disp(['5 shells: Vw = ' num2str(vw(2))]);
%
vw(1) = rep_4sh;
vw(2) = rep_5sh;



% Plot
bar(1:n_ep, vw, 'facecolor', bar_color, 'linewidth', lw_bar, 'barwidth', bw1)

% Axis and labels
set(gca, 'xtick', 1:n_ep, 'xticklabel', ep_list_shell_name, ...
    'box', 'off', 'tickdir', 'out', 'fontweight', 'bold', 'linewidth', ...
    lw_axis, 'fontsize', fs_axis, 'ticklength', tl_axis, ...
     'ytick', ytick, 'yticklabel', yticklbl)
xlim([0 n_ep+1])
ylim([0 ymax_vw])
xlabel(shell_lbl, 'fontsize', fs_lbl)
ylabel(sdw_str, 'fontsize', fs_lbl)
%
yl = get(gca, 'ylabel');
p  = get(yl, 'position');
set(yl, 'position', [p(1) yl_pos p(3)])


%%% PRINT FIGURE TO FILE
output_filename = fullfile(output_dir, output_name);
save_current_fig_to_file(output_name, output_dir, imsize, res_dpi);
disp(['Wrote: ' output_filename]);
system(['open ' output_filename]);

end






