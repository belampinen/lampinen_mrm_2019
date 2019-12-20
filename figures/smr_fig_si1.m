function smr_fig_si1(output_dir)

if (nargin < 1), output_dir = fullfile(pwd, 'figures', 'output'); end

output_name = 'smr_fig_si1.tiff';

%%%------------------------------------------------------------------
% PREPARE
%%%------------------------------------------------------------------

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
n_ep = numel(ep_list);
%
ep_name_bdelta = {...
    '1', ...
    '1/-0.5', ...
    '1/0', ...
    'any'};
%
ep_name_shell = {...    
    '4', ...
    '5', ...
    '6', ...
    '7', ...
    '8', ...
    '9', ...
    '10'};
%
opt_optimize        = smr_optimize_opt_derive(smr_optimize_opt);
opt_optimize.par_sign_diff(5) = 0.1;
%
%
c_ep_bdelta = [1 2 3 10];
n_ep_bdelta = numel(c_ep_bdelta);
c_ep_shell  = 4:10;
n_ep_shell  = numel(c_ep_shell);

%
yticklbl        = {'0.0', '1.0', '2.0', '3.0', '4.0', '5.0'};
tl_axis = 0.032 + [0 0];

%
bootstrap_name = 'Simulation';
%
a_str           = 'Shapes of the b-tensor';
b_str           = 'Number of shells';
%
vw_str         = 'V_{W}';
%
shell_lbl           = 'Number of shells';
b_delta_lbl = '\bf\itb\rm\bf_{\Delta} combinations';


% Plotting
%
n_cols_figure   = 1;
journal         = 'mrm';
res_dpi         = 600;
%
%
%
m_left   = 0.35;
pw       = 1;
pwm      = 0.6;
m_right  = 0.25;
%
m_upper = 0.05;
th      = 0.2;
thm     = 0.25;  % for broken bars
ph      = 1;
m_lower = 0.4;
%
tm      = 0.25;
tw      = 2; % A, B, C, D...
%
axis_factor       = 1;
leg_start_part_l  = 0.7;
leg_start_part_b  = 0.6;

%

% Axis setup
fh = m_upper    + th + thm + ph + m_lower;
fw = m_left     + pw + pwm + pw + m_right;
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
fs_title        = 7;
fs_other        = 5;
fs_legend       = fs_other;
fs_axis         = fs_other;
fs_lbl          = fs_other;
%
lw_axis         = 1;
lw_bar          = 0.5;
%
bar_color       = 0.8 + [0 0 0];
bw1             = 0.3;
bw2             = 0.2;
bw_scale_leg    = 0.3;
%
ms_crlb_k       = 8;
%
ymax_vw         = 5.5;
%
% CRLB-plots
lw_crlb          = lw_bar;
its_crlb         = 12;


%%%------------------------------------------------------------------
% PERFORM
%%%------------------------------------------------------------------

%%% Obtain vw values
%
% -- Vw from CRLB
vw_crlb = zeros(n_ep, 1);
for c_ep = 1:n_ep
    vw_crlb(c_ep) = smr_optimize_ep2vw(smr_ep(ep_list{c_ep}), opt_optimize);
end
%
% -- Vw from Simulation
mat_fn    = fullfile(output_dir, 'smr_fig_si1_sim.mat');
my_struct = load(mat_fn);  
vw_sim    = my_struct.vw;


disp(['Vw CRLB for LTE-only = ' num2str(vw_crlb(1))]);
disp(['Vw CRLB for 4 shells = ' num2str(vw_crlb(4))]);
disp(['Vw CRLB for 5 shells = ' num2str(vw_crlb(5))]);
%
disp(['Vw Sim for LTE-only = ' num2str(vw_sim(1))]);
disp(['Vw Sim for 4 shells = ' num2str(vw_sim(4))]);
disp(['Vw Sim for 5 shells = ' num2str(vw_sim(5))]);

%%% Figure
%
figure(967)
clf
set(gcf, 'color', 'w')


%%% ------- A: Shapes of the B-tensor ------- %%%
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

%

% Legend
bar (-1, 0, 'facecolor', bar_color, 'linewidth', lw_bar, 'barwidth', bw_scale_leg * bw2);
hold on;
plot(-1,0, 'pk', 'markersize', ms_crlb_k, 'linewidth', lw_crlb)
l = legend({'CRLB', bootstrap_name}, 'box', 'off', 'fontsize', fs_legend);
pos = get(l, 'pos');
set(l, 'pos', [(ax_l + leg_start_part_l * pw) (ax_b + leg_start_part_b * ph) pos(3:4)]) % pos(3) pos(4)]) % lbwh
set(l, 'ItemTokenSize', [10 its_crlb]);

% Plot
bar( 1:n_ep_bdelta, vw_crlb(c_ep_bdelta), 'facecolor', bar_color, 'linewidth', lw_bar, 'barwidth', bw2)
hold on
plot(1:n_ep_bdelta, vw_sim(c_ep_bdelta) , 'pk', 'markersize', ms_crlb_k, 'linewidth', lw_crlb)

% Axis and labels
set(gca, 'xtick', 1:n_ep_bdelta, 'xticklabel', ep_name_bdelta, ...
    'box', 'off', 'tickdir', 'out', 'fontweight', 'bold', 'linewidth', lw_axis, 'fontsize', fs_axis, ...
    'xticklabelrot', 0,  'ticklength', tl_axis, 'yticklabel', yticklbl)
xlim([0 n_ep_bdelta+1])
ylim([0 ymax_vw])
if (~isempty(b_delta_lbl))
    xlabel(b_delta_lbl, 'fontsize', fs_lbl)
end
ylabel(vw_str, 'fontsize', fs_lbl)
axis square



%%% ------- B: Combinations ------- %%%
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


% Legend
bar (-1, 0, 'facecolor', bar_color, 'linewidth', lw_bar, 'barwidth', bw_scale_leg * bw1);
hold on;
plot(-1,0, 'pk', 'markersize', ms_crlb_k, 'linewidth', lw_crlb)
l = legend({'CRLB', bootstrap_name}, 'box', 'off', 'fontsize', fs_legend);
pos = get(l, 'pos');
set(l, 'pos', [(ax_l + leg_start_part_l * pw) (ax_b + leg_start_part_b * ph) pos(3:4)]) % pos(3) pos(4)]) % lbwh
set(l, 'ItemTokenSize', [10 its_crlb]);

% Plot
bar (1:n_ep_shell, vw_crlb(c_ep_shell), 'facecolor', bar_color, 'linewidth', lw_bar, 'barwidth', bw1)
hold on
plot(1:n_ep_shell, vw_sim(c_ep_shell), 'pk', 'markersize', ms_crlb_k, 'linewidth', lw_crlb)

% Axis and labels
set(gca, 'xtick', 1:n_ep_shell, 'xticklabel', ep_name_shell, ...
    'box', 'off', 'tickdir', 'out', 'fontweight', 'bold', 'linewidth', lw_axis, 'fontsize', fs_axis, ...
     'ticklength', tl_axis, 'yticklabel', yticklbl)
xlim([0 n_ep_shell+1])
ylim([0 ymax_vw])
xlabel(shell_lbl, 'fontsize', fs_lbl)
ylabel(vw_str, 'fontsize', fs_lbl)
axis square



%%% PRINT FIGURE TO FILE
output_filename = fullfile(output_dir, output_name);
save_current_fig_to_file(output_name, output_dir, imsize, res_dpi);
disp(['Wrote: ' output_filename]);
system(['open ' output_filename]);

end