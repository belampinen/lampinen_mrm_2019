function smr_fig6(output_dir)

if (nargin < 1), output_dir = fullfile(pwd, 'figures', 'output'); end

output_name             = 'smr_fig6.tiff';


%%%------------------------------------------------------------------
% PREPARE
%%%------------------------------------------------------------------

do_legend = 0;
%
c_roi_wm = [2 3 4 5];
c_roi_gm = [6 7 8 9];
%
% roi_name_list = {...
%     'noncortex', ...
%     ...
%     'cso', ...
%     'plic', ...
%     'or', ...
%     'acr', ...
%     ...
%     'thal', ...
%     'caud', ...
%     'put', ...
%     'gp', ...
%     };
%
par_name_list = {...
    [], ...
    '\bf\itf\rm\bf_{S}', ...                  % 2  +
    '\bf\itD\rm\bf_{I;S} [\mum^2/ms]', ...    % 3
    '\bf\itD\rm\bf_{I;Z} [\mum^2/ms]', ...    % 4
    '\bf\itD\rm\bf_{\Delta;Z}', ...           % 5  +
    [], ...
    [], ...
    [], ...
    [], ...
    []', ...
    '\bfT_{2;S} [ms]', ...                  % 11
    '\bfT_{2;Z} [ms]', ...                  % 12
    '\bf\itp\rm\bf_{2}', ...               % 13
    [], ...                %
    '\bf\itD\rm\bf_{||;Z} [\mum^2/ms]', ...   % 15
    '\bf\itD\rm\bf_{\perp;Z} [\mum^2/ms]', ...% 16
    [], ...
    '\bf\itD\rm\bf_{||;S} [\mum^2/ms]', ...   % 18 
    };

t2_max = 200;
par_range_list = {...
    [], [0 1], [0 1], [0 2], [0 1], [], [], [], [], [], [0 t2_max], [0 t2_max], [0 1], [], ...
    [0 4], [0 4], [], [0 4]};
%
diff_tick_lbl = {'0', '1', '2', '3'}; % , '4'};
ticklbl_list = {...
    [], ...
    {'0.00', '0.25', '0.50', '0.75', '1.00'}, ...
    {'0.00', '0.25', '0.50', '0.75', '1.00'}, ...
    {'0.0', '0.5', '1.0', '1.5', '2.0'}, ...
    {'0.00', '0.25', '0.50', '0.75', '1.00',}, ...
    [], ...
    [], ...
    [], ...
    [], ...
    [], ...
    {'0', '50', '100', '150', '200'}, ...
    {'0', '50', '100', '150', '200'}, ...
    {'0.00', '0.25', '0.50', '0.75', '1.00'}, ...
    [], ...
    diff_tick_lbl, ...
    diff_tick_lbl, ...
    [], ...
    diff_tick_lbl, ...
    };

%%%
%
roi_col_wm   = [0 0 0];
roi_col_gm   = [0.3 0.6 0.3];
roi_col_wml  = [0.9 0.2 0.2];
%
dens_txt_x = 135;
dens_txt_y = 160;
%
propdiff_txt_x1 = 0.35;
propdiff_txt_y1 = 2.3;
propdiff_txt_x2 = 1.4;
propdiff_txt_y2 = 1.95;
%
tort_txt_x = 0.42;
tort_txt_y = 0.06;
tort_rot   = 35;
%
tl_axis = 0.025 + [0 0];


%
a_str           = 'A - Parameter relations commonly assumed in modeling';
b_str           = 'B - Parameter relations potentially supported by data';
%
% tp_sc = 1.3;
% tpy = 3.9625;
%
a1_str = 'Tortuosity constraint';
a2_str = 'Proportional diffusivities';
a3_str = 'Density assumption';
%
b1_str = 'Shape \propto coherence';
b2_str = 'Diffusivity \propto fraction';
b3_str = 'Tortuosity with g-ratio';
%
plot_order_a = {'a3', 'a2', 'a1'};

%%%------------------------------------------------------------------
% PERFORM
%%%------------------------------------------------------------------

%%% Load values
%

my_struct = load('estimates.mat');

%             s0  fs dis diz ddz ---SH---   t2s t2z  msr p2 daz drz t2  das g mwf
sc          = [1  1  1   1   1   1 1 1 1 1  1   1    1   1  1   1   1   1   1 1];
%
vals_wmgm = squeeze(my_struct.fit_vals_m(1:10,:,:));
vals_wml = squeeze(my_struct.fit_vals_m(11:15,5,:));

    


%%% Plot
%
figure_nr       = 278;
n_cols_figure   = 2;
journal         = 'mrm';
res_dpi         = 600;
%
m_upper = 0.07;
th      = 0.2;
thm     = 0.0;
th2     = 0.15;
thm2    = 0.05;
ph      = 1;
iph     = 0.45;
m_lower = 0.35;
%
m_left   = 0.38;
pw       = 1;
pwm      = 0.5;
m_right  = 0.2;
%
tm      = 0.32;
tw      = 2; % A and B
%
axis_factor     = 1;
%
n_plot = 3;

% Axis setup
fh = m_upper + 2 * (th + thm + th2 + thm2 + ph) + iph + m_lower;
fw = m_left  + n_plot * pw + (n_plot - 1) * pwm + m_right;
%
m_upper = m_upper / fh;
th      = th / fh;
thm     = thm / fh;
th2     = th2 / fh;
thm2    = thm2 / fh;
ph      = ph / fh;
iph     = iph / fh;
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
%
my_fw       = 'bold';
%
fs_txt          = 8;
fs_axis         = 9;
fs_label        = 9;
fs_title_panel  = 11;
fs_title        = 9;
%
lw_axis     = 1;
lw_con      = 1;
lw_g        = 0.6;
%
n_step = 20;
%
dot_size    = 18;
fa          = 0.4;
lw_scatter  = 0.7;

%%% Plot
%
figure(figure_nr)
clf
set(gcf, 'color', 'w')

%


%%% ------- Panel A ------- %%%
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
    'fontweight', my_fw, ...
    'FitBoxToText', 'on', ...
    'linestyle', 'none', ...
    'horizontalalignment', 'left', ...
    'verticalalignment', 'middle', ...
    'fontsize', fs_title_panel);
%
%
for c_plot = 1:numel(plot_order_a)
    
    switch (plot_order_a{c_plot})
        
        case 'a1'
            
            %%% 1: Tortuosity
            %
            %
            % A1 annotation
            t_l = l_pos;
            t_b = 1 - m_upper - th - thm - th2;
            t_w = pw;
            t_h = th2;
            annotation(...
                'textbox', [t_l t_b t_w t_h], ...
                'string', a1_str, ...
                'fontweight', my_fw, ...
                'FitBoxToText', 'on', ...
                'linestyle', 'none', ...
                'horizontalalignment', 'center', ...
                'verticalalignment', 'middle', ...
                'fontsize', fs_title);
            %
            %
            ax_l = l_pos;
            ax_b = 1 - m_upper - th - thm - th2 - thm2 - ph;
            ax_w = pw * axis_factor;
            ax_h = ph * axis_factor;
            af_adjust = - pw * (axis_factor - 1) / 2; % Correct for axis factor
            ax_l = ax_l + af_adjust;
            ax_b = ax_b + af_adjust;
            axes('position', [ax_l ax_b ax_w ax_h]);
            %
            l_pos = l_pos + pw + pwm;
            
            %
            c_par_x = 2;  % fs
            c_par_y = 5;  % dd_z
            %
            
            % Plot constraint
            %
            %
            xl = par_range_list{c_par_x};
            yl = par_range_list{c_par_y};
            x  = linspace(xl(1), xl(2), n_step);
            y  = x ./ (3 - 2 * x);
            %
            plot(x,y,'--k', 'linewidth', lw_con)
            hold on
            %
            txt = '\bf\itD\rm\bf_{\Delta;Z} = \itf\rm\bf_{S}/(3 - 2\itf\rm\bf_{S})';
            %
            
            text(tort_txt_x, tort_txt_y, txt, 'color', 'k', ...
                'fontsize', fs_txt, 'rotation', tort_rot)
            
            
            % Extract values
            1;
            
            [x_wm, x_gm, x_wml, y_wm, y_gm, y_wml] = ...
                get_vals(vals_wmgm, vals_wml, c_roi_wm, c_roi_gm, c_par_x, c_par_y, sc);
            
            % Plot values
            scatter(x_gm, y_gm, dot_size, ...
                'markeredgecolor', roi_col_gm, 'markerfacecolor', roi_col_gm, ...
                'MarkerFaceAlpha', fa, 'linewidth', lw_scatter);
            hold on 
            scatter(x_wm, y_wm, dot_size, ...
                'markeredgecolor', roi_col_wm, 'markerfacecolor', roi_col_wm, ...
                'MarkerFaceAlpha', fa, 'linewidth', lw_scatter);
            scatter(x_wml, y_wml, dot_size, ...
                'markeredgecolor', roi_col_wml, 'markerfacecolor', roi_col_wm, ...
                'MarkerFaceAlpha', fa, 'linewidth', lw_scatter);
            
            axis square;
            
            %
            set(gca, 'tickdir', 'out', 'fontsize', fs_axis, 'box', 'off', ...
                'linewidth', lw_axis, 'fontweight', my_fw, 'ticklength', tl_axis)
            %
            xlim(xl);
            ylim(yl);
            %
            set(gca, 'xtick',  linspace(xl(1), xl(2), numel(ticklbl_list{c_par_x})), 'xticklabel', ticklbl_list{c_par_x});
            set(gca, 'ytick',  linspace(yl(1), yl(2), numel(ticklbl_list{c_par_y})), 'yticklabel', ticklbl_list{c_par_y});
            %
            xlabel(par_name_list{c_par_x}, 'fontsize', fs_label)
            ylabel(par_name_list{c_par_y}, 'fontsize', fs_label)
            
            
        case 'a2'
            
            %%% 2: Diffusivity
            %
            % A2 annotation
            t_l = l_pos;
            t_b = 1 - m_upper - th - thm - th2;
            t_w = pw;
            t_h = th2;
            annotation(...
                'textbox', [t_l t_b t_w t_h], ...
                'string', a2_str, ...
                'fontweight', my_fw, ...
                'FitBoxToText', 'on', ...
                'linestyle', 'none', ...
                'horizontalalignment', 'center', ...
                'verticalalignment', 'middle', ...
                'fontsize', fs_title);
            %
            %
            ax_l = l_pos;
            ax_b = 1 - m_upper - th - thm - th2 - thm2 - ph;
            ax_w = pw * axis_factor;
            ax_h = ph * axis_factor;
            af_adjust = - pw * (axis_factor - 1) / 2; % Correct for axis factor
            ax_l = ax_l + af_adjust;
            ax_b = ax_b + af_adjust;
            axes('position', [ax_l ax_b ax_w ax_h]);
            %
            l_pos = l_pos + pw + pwm;
            
            %
            c_par_x = 15;  % da_z
            c_par_y = 18;  % da_s
            %
            
            % Plot constraint
            %
            %
            xl = par_range_list{c_par_x};
            yl = par_range_list{c_par_y};
            x  = linspace(xl(1), xl(2), n_step);
            %
            y1 = 1 * x; % NODDI/SMT/B&S
            %
            plot(x , y1, '--k', 'linewidth', lw_con)
            hold on
            %                               
            txt = '\bf\itD\rm\bf_{||;S} = \itD\rm\bf_{||;Z}';
            text(propdiff_txt_x2, propdiff_txt_y2, txt, 'color', 'k', ...
                'fontsize', fs_txt, 'rotation', 45);
            %
            y2 = 3 * x; % CODIVIDE/rNODDI
            %
            plot(x , y2, '--k', 'linewidth', lw_con)
            %
            txt = '\bf\itD\rm\bf_{||;S} = 3\itD\rm\bf_{||;Z}';
            text(propdiff_txt_x1, propdiff_txt_y1, txt, 'color', 'k', ...
                'fontsize', fs_txt, 'rotation', 72);
            %
            
            % Extract values
            [x_wm, x_gm, x_wml, y_wm, y_gm, y_wml] = ...
                get_vals(vals_wmgm, vals_wml, c_roi_wm, c_roi_gm, c_par_x, c_par_y, sc);
            
            
            % Plot values
            scatter(x_wm, y_wm, dot_size, ...
                'markeredgecolor', roi_col_wm, 'markerfacecolor', roi_col_wm, ...
                'MarkerFaceAlpha', fa, 'linewidth', lw_scatter);
            hold on
            scatter(x_gm, y_gm, dot_size, ...
                'markeredgecolor', roi_col_gm, 'markerfacecolor', roi_col_gm, ...
                'MarkerFaceAlpha', fa, 'linewidth', lw_scatter);
            %
            scatter(x_wml, y_wml, dot_size, ...
                'markeredgecolor', roi_col_wml, 'markerfacecolor', roi_col_wm, ...
                'MarkerFaceAlpha', fa, 'linewidth', lw_scatter);
            
            %
            axis square;
            %
            set(gca, 'tickdir', 'out', 'fontsize', fs_axis, 'box', 'off', ...
                'linewidth', lw_axis, 'fontweight', my_fw, 'ticklength', tl_axis)
            %
            xlim(xl);
            ylim(yl);
            %
            set(gca, 'xtick',  linspace(xl(1), xl(2), numel(ticklbl_list{c_par_x})), 'xticklabel', ticklbl_list{c_par_x});
            set(gca, 'ytick',  linspace(yl(1), yl(2), numel(ticklbl_list{c_par_y})), 'yticklabel', ticklbl_list{c_par_y});
            %
            xlabel(par_name_list{c_par_x}, 'fontsize', fs_label)
            ylabel(par_name_list{c_par_y}, 'fontsize', fs_label)
       
            
        case 'a3'
            
            %%% 3: Density
            %
            % A3 annotation
            t_l = l_pos;
            t_b = 1 - m_upper - th - thm - th2;
            t_w = pw;
            t_h = th2;
            annotation(...
                'textbox', [t_l t_b t_w t_h], ...
                'string', a3_str, ...
                'fontweight', my_fw, ...
                'FitBoxToText', 'on', ...
                'linestyle', 'none', ...
                'horizontalalignment', 'center', ...
                'verticalalignment', 'middle', ...
                'fontsize', fs_title);
            %
            %
            ax_l = l_pos;
            ax_b = 1 - m_upper - th - thm - th2 - thm2 - ph;
            ax_w = pw * axis_factor;
            ax_h = ph * axis_factor;
            af_adjust = - pw * (axis_factor - 1) / 2; % Correct for axis factor
            ax_l = ax_l + af_adjust;
            ax_b = ax_b + af_adjust;
            axes('position', [ax_l ax_b ax_w ax_h]);
            %
            l_pos = l_pos + pw + pwm;
            
            %
            c_par_x = 11;  % T2s
            c_par_y = 12;  % T2z
            %
            
            % Plot constraint
            %
            %
            xl = par_range_list{c_par_x};
            yl = par_range_list{c_par_y};
            x  = linspace(xl(1), xl(2), n_step);
            %
            y  = x; % Density
            %
            plot(x , y, '--k', 'linewidth', lw_con)
            hold on
            %
            txt = '\bfT2_Z = T2_S';
            text(dens_txt_x, dens_txt_y, txt, 'color', 'k', ...
                'fontsize', fs_txt, 'rotation', 45);
            %           
            
            % Extract values
            [x_wm, x_gm, x_wml, y_wm, y_gm, y_wml] = ...
                get_vals(vals_wmgm, vals_wml, c_roi_wm, c_roi_gm, c_par_x, c_par_y, sc);
            
            % Plot values
            scatter(x_gm, y_gm, dot_size, ...
                'markeredgecolor', roi_col_gm, 'markerfacecolor', roi_col_gm, ...
                'MarkerFaceAlpha', fa, 'linewidth', lw_scatter);
            hold on
            scatter(x_wm, y_wm, dot_size, ...
                'markeredgecolor', roi_col_wm, 'markerfacecolor', roi_col_wm, ...
                'MarkerFaceAlpha', fa, 'linewidth', lw_scatter);
            %
            scatter(x_wml, y_wml, dot_size, ...
                'markeredgecolor', roi_col_wml, 'markerfacecolor', roi_col_wm, ...
                'MarkerFaceAlpha', fa, 'linewidth', lw_scatter);
            
            %
            axis square;
            %
            set(gca, 'tickdir', 'out', 'fontsize', fs_axis, 'box', 'off', ...
                'linewidth', lw_axis, 'fontweight', my_fw, 'ticklength', tl_axis)
            %
            xlim(xl);
            ylim(yl);
            %
            set(gca, 'xtick',  linspace(xl(1), xl(2), numel(ticklbl_list{c_par_x})), 'xticklabel', ticklbl_list{c_par_x});
            set(gca, 'ytick',  linspace(yl(1), yl(2), numel(ticklbl_list{c_par_y})), 'yticklabel', ticklbl_list{c_par_y});
            %
            xlabel(par_name_list{c_par_x}, 'fontsize', fs_label)
            ylabel(par_name_list{c_par_y}, 'fontsize', fs_label)
            
            
    end
end


%%% ------- Panel B ------- %%%
%
l_pos = m_left;
row_a_h = th + thm + th2 + thm2 + ph + iph;
%
% B annotation
t_l = l_pos - tm;
t_b = 1 - m_upper - th - row_a_h;
t_w = pw;
t_h = th;
annotation(...
    'textbox', [t_l t_b t_w t_h], ...
    'string', b_str, ...
    'fontweight', my_fw, ...
    'FitBoxToText', 'on', ...
    'linestyle', 'none', ...
    'horizontalalignment', 'left', ...
    'verticalalignment', 'middle', ...
    'fontsize', fs_title_panel);
%


%%% 1: Ddelta / p2
%
% B1 annotation
t_l = l_pos;
t_b = 1 - m_upper - th - thm - th2 - row_a_h;
t_w = pw;
t_h = th2;
annotation(...
    'textbox', [t_l t_b t_w t_h], ...
    'string', b1_str, ...
    'fontweight', my_fw, ...
    'FitBoxToText', 'on', ...
    'linestyle', 'none', ...
    'horizontalalignment', 'center', ...
    'verticalalignment', 'middle', ...
    'fontsize', fs_title);
%
%
ax_l = l_pos;
ax_b = 1 - m_upper - th - thm - th2 - thm2 - ph - row_a_h;
ax_w = pw * axis_factor;
ax_h = ph * axis_factor;
af_adjust = - pw * (axis_factor - 1) / 2; % Correct for axis factor
ax_l = ax_l + af_adjust;
ax_b = ax_b + af_adjust;
axes('position', [ax_l ax_b ax_w ax_h]);
%
l_pos = l_pos + pw + pwm;


%
c_par_x = 13;         % p2
c_par_y = 5;          % dd_z
%
% Extract values
[x_wm, x_gm, x_wml, y_wm, y_gm, y_wml] = ...
    get_vals(vals_wmgm, vals_wml, c_roi_wm, c_roi_gm, c_par_x, c_par_y, sc);

% Plot values
scatter(x_gm, y_gm, dot_size, ...
    'markeredgecolor', roi_col_gm, 'markerfacecolor', roi_col_gm, ...
    'MarkerFaceAlpha', fa, 'linewidth', lw_scatter);
hold on
scatter(x_wm, y_wm, dot_size, ...
    'markeredgecolor', roi_col_wm, 'markerfacecolor', roi_col_wm, ...
    'MarkerFaceAlpha', fa, 'linewidth', lw_scatter);
%
scatter(x_wml, y_wml, dot_size, ...
    'markeredgecolor', roi_col_wml, 'markerfacecolor', roi_col_wm, ...
    'MarkerFaceAlpha', fa, 'linewidth', lw_scatter);
%
axis square;
%
set(gca, 'tickdir', 'out', 'fontsize', fs_axis, 'box', 'off', ...
    'linewidth', lw_axis, 'fontweight', my_fw, 'ticklength', tl_axis)
%
xl = par_range_list{13};
yl = par_range_list{c_par_y};
xlim(xl);
ylim(yl);
%
set(gca, 'xtick',  linspace(xl(1), xl(2), numel(ticklbl_list{c_par_x})), 'xticklabel', ticklbl_list{c_par_x});
set(gca, 'ytick',  linspace(yl(1), yl(2), numel(ticklbl_list{c_par_y})), 'yticklabel', ticklbl_list{c_par_y});
%
xlabel(par_name_list{c_par_x}, 'fontsize', fs_label)
ylabel(par_name_list{c_par_y}, 'fontsize', fs_label)


%%% 2: Dis/fs
%
% B2 annotation
t_l = l_pos;
t_b = 1 - m_upper - th - thm - th2 - row_a_h;
t_w = pw;
t_h = th2;
annotation(...
    'textbox', [t_l t_b t_w t_h], ...
    'string', b2_str, ...
    'fontweight', my_fw, ...
    'FitBoxToText', 'on', ...
    'linestyle', 'none', ...
    'horizontalalignment', 'center', ...
    'verticalalignment', 'middle', ...
    'fontsize', fs_title);
%
%
ax_l = l_pos;
ax_b = 1 - m_upper - th - thm - th2 - thm2 - ph - row_a_h;
ax_w = pw * axis_factor;
ax_h = ph * axis_factor;
af_adjust = - pw * (axis_factor - 1) / 2; % Correct for axis factor
ax_l = ax_l + af_adjust;
ax_b = ax_b + af_adjust;
axes('position', [ax_l ax_b ax_w ax_h]);
%
l_pos = l_pos + pw + pwm;


%
c_par_x = 2;    % fs
c_par_y = 18;   % da_s
%
% Extract values
[x_wm, x_gm, x_wml, y_wm, y_gm, y_wml] = ...
    get_vals(vals_wmgm, vals_wml, c_roi_wm, c_roi_gm, c_par_x, c_par_y, sc);

% Plot values
scatter(x_gm, y_gm, dot_size, ...
    'markeredgecolor', roi_col_gm, 'markerfacecolor', roi_col_gm, ...
    'MarkerFaceAlpha', fa, 'linewidth', lw_scatter);
hold on
scatter(x_wm, y_wm, dot_size, ...
    'markeredgecolor', roi_col_wm, 'markerfacecolor', roi_col_wm, ...
    'MarkerFaceAlpha', fa, 'linewidth', lw_scatter);
%
scatter(x_wml, y_wml, dot_size, ...
    'markeredgecolor', roi_col_wml, 'markerfacecolor', roi_col_wm, ...
    'MarkerFaceAlpha', fa, 'linewidth', lw_scatter);
%
axis square;
%
set(gca, 'tickdir', 'out', 'fontsize', fs_axis, 'box', 'off', ...
    'linewidth', lw_axis, 'fontweight', my_fw, 'ticklength', tl_axis)
%
xl = par_range_list{c_par_x};
yl = par_range_list{c_par_y};
xlim(xl);
ylim(yl);
%
set(gca, 'xtick',  linspace(xl(1), xl(2), numel(ticklbl_list{c_par_x})), 'xticklabel', ticklbl_list{c_par_x});
set(gca, 'ytick',  linspace(yl(1), yl(2), numel(ticklbl_list{c_par_y})), 'yticklabel', ticklbl_list{c_par_y});
%
xlabel(par_name_list{c_par_x}, 'fontsize', fs_label)
ylabel(par_name_list{c_par_y}, 'fontsize', fs_label)


% B3 annotation
t_l = l_pos;
t_b = 1 - m_upper - th - thm - th2 - row_a_h;
t_w = pw;
t_h = th2;
annotation(...
    'textbox', [t_l t_b t_w t_h], ...
    'string', b3_str, ...
    'fontweight', my_fw, ...
    'FitBoxToText', 'on', ...
    'linestyle', 'none', ...
    'horizontalalignment', 'center', ...
    'verticalalignment', 'middle', ...
    'fontsize', fs_title);
%

col_list = {...
    [41 73 203 ]/255, ...
    [41 203 203]/255, ...
    [144 203 41]/255, ...
    [203 171 41]/255, ...
    [203 84 41 ]/255, ...
    }; %
%
%
ax_l = l_pos;
ax_b = 1 - m_upper - th - thm - th2 - thm2 - ph - row_a_h;
ax_w = pw * axis_factor;
ax_h = ph * axis_factor;
af_adjust = - pw * (axis_factor - 1) / 2; % Correct for axis factor
ax_l = ax_l + af_adjust;
ax_b = ax_b + af_adjust;
axes('position', [ax_l ax_b ax_w ax_h]);
%
%
c_par_x = 2;  % fs
c_par_y = 5;  % dd_z
%

% Legend
for c_g = 1:numel(col_list)
    plot(-1, -1, '--', 'linewidth', lw_g, 'color', col_list{c_g});
    hold on
end

%
xl = par_range_list{c_par_x};
yl = par_range_list{c_par_y};

%%% Plot 'corrected' constraints using different g-ratios
%
g_vec   = linspace(0.4, 0.99, 5);
% [165 41 203]/255};
x       = linspace(0, 1, n_step);
f_fsp   = @(g) x * (1 + (g^2 / (1 - g^2))) ./ (x + (g^2 / (1 - g^2)));
%
for c_g = 1:numel(g_vec)
    
    fsp = f_fsp(g_vec(c_g));
    y   = fsp ./ (3 - 2 * fsp);
    
    plot(x, y, '--', 'color', col_list{c_g}, 'linewidth', lw_g)
    hold on
end
%
set(gca, 'tickdir', 'out', 'fontsize', fs_axis, 'box', 'off', ...
    'linewidth', lw_axis, 'fontweight', my_fw, 'ticklength', tl_axis)
%
%
txt = '\bf\itD\rm\bf_{\Delta;Z} = \itv\rm\bf_{E}/(3 - 2\itv\rm\bf_{E})';

text(tort_txt_x, tort_txt_y, txt, 'color', 'k', ...
    'fontsize', fs_txt, 'rotation', tort_rot);
%
% Extract values
[x_wm, x_gm, x_wml, y_wm, y_gm, y_wml] = ...
    get_vals(vals_wmgm, vals_wml, c_roi_wm, c_roi_gm, c_par_x, c_par_y, sc);

% Plot values
scatter(x_wm, y_wm, dot_size, ...
    'markeredgecolor', roi_col_wm, 'markerfacecolor', roi_col_wm, ...
    'MarkerFaceAlpha', fa, 'linewidth', lw_scatter);
hold on
scatter(x_wml, y_wml, dot_size, ...
    'markeredgecolor', roi_col_wml, 'markerfacecolor', roi_col_wm, ...
    'MarkerFaceAlpha', fa, 'linewidth', lw_scatter);
axis square;


xlim(xl);
ylim(yl);
%
set(gca, 'xtick',  linspace(xl(1), xl(2), numel(ticklbl_list{c_par_x})), 'xticklabel', ticklbl_list{c_par_x});
set(gca, 'ytick',  linspace(yl(1), yl(2), numel(ticklbl_list{c_par_y})), 'yticklabel', ticklbl_list{c_par_y});
%
xlabel(par_name_list{c_par_x}, 'fontsize', fs_label)
ylabel(par_name_list{c_par_y}, 'fontsize', fs_label)


if (do_legend)
    fs_leg = 5;
    leg_str = {'\bf\itg\rm\bf = 0.4', '\bf\itg\rm\bf = 0.55', '\bf\itg\rm\bf = 0.7', '\bf\itg\rm\bf = 0.85', '\bf\itg\rm\bf = 1'};
    l = legend(leg_str, 'box', 'off', 'fontsize', fs_leg);
    p = get(l, 'pos');
    set(l, 'pos', [1.104 * p(1) 2 * p(2) p(3) p(4)])
end

%%% PRINT FIGURE TO FILE
output_filename = fullfile(output_dir, output_name);
save_current_fig_to_file(output_name, output_dir, imsize, res_dpi);
disp(['Wrote: ' output_filename]);
system(['open ' output_filename]);

end


%------------------------------------------------------------------------------------
function [x_wm, x_gm, x_wml, y_wm, y_gm, y_wml] ...
    = get_vals(vals_wmgm, vals_wml, c_roi_wm, c_roi_gm, c_par_x, c_par_y, sc)

tmp        = squeeze(vals_wmgm(:,c_roi_wm,c_par_x));
x_wm       = tmp(:) * sc(c_par_x);
tmp        = squeeze(vals_wmgm(:,c_roi_gm,c_par_x));
x_gm       = tmp(:) * sc(c_par_x);
x_wml      = vals_wml(:,c_par_x) * sc(c_par_x);

%
tmp        = squeeze(vals_wmgm(:,c_roi_wm,c_par_y));
y_wm       = tmp(:) * sc(c_par_y);
tmp        = squeeze(vals_wmgm(:,c_roi_gm,c_par_y));
y_gm       = tmp(:) * sc(c_par_y);
y_wml      = vals_wml(:,c_par_y) * sc(c_par_y);

end
