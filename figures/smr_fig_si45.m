function smr_fig_si45(output_dir)

if (nargin < 1), output_dir = fullfile(pwd, 'figures', 'output'); end


if (0)
    smr_fig_si45_sim;
end


%%%------------------------------------------------------------------
% PREPARE
%%%------------------------------------------------------------------

mp_name_full_od = 'a_full_od';
mp_name_mid_od  = 'a_mid_od';
mp_list = {...
    mp_name_full_od, ...
    mp_name_full_od, ...
    mp_name_full_od, ...
    mp_name_full_od, ...
    mp_name_mid_od, ...
    mp_name_mid_od, ...
    mp_name_mid_od, ...
    mp_name_mid_od};
%
my_green    = [0.40 0.70 0.40];
my_yellow   = [0.98 0.83 0.00];
my_red      = [0.70 0.30 0.30];
%
ep_col_list      = {...
    my_green, ...
    my_yellow, ...
    my_red, ...
    'k', ...
    my_green, ...
    my_yellow, ...
    my_red, ...
    'k', ...
    };
    
% Choose plot case
if (0) 
    c_case_plot     = [4 3 2 1];
    output_name     = 'smr_fig_si4.tiff';
    gt_p2           = 0;
else  
    c_case_plot     = [8 7 6 5];
    output_name     = 'smr_fig_si5.tiff';
    gt_p2           = 0.14;
end


%
par_name_list = {...
    's0', ...
    '\bf\itf\rm\bf_{S}', ...
    ['\bf\itD\rm\bf_{I;S} [' char(181) 'm^2/ms]'], ...
    ['\bf\itD\rm\bf_{I;Z} [' char(181) 'm^2/ms]'], ...
    '\bf\itD\rm\bf_{\Delta;Z}', ...
    'p20', ...
    'p21r', ...
    'p21i', ...
    'p22r', ...
    'p22i', ...
    '\bfT_{2;S} [ms]', ...
    '\bfT_{2;Z} [ms]', ...    
    '\bf\itp\rm\bf_{2}', ...  
    };

%
t2_max = 150;
t2_min = 0;
di_max = 3;
par_range_list = {...
    [], [0 1], [0 di_max], [0  di_max], [-0.5 1], [], [], [], [], [], [t2_min t2_max], [t2_min t2_max], [0 1]};
%
ticklbl_list = {...
    [], ...
    {'0', '0.5', '1'}, ...
    {'0', num2str(di_max/2), num2str(di_max)}, ...
    {'0', num2str(di_max/2), num2str(di_max)}, ...
    {'-0.5', '0', '0.5', '1'}, ...
    [], ...
    [], ...
    [], ...
    [], ...
    [], ...
    {num2str(t2_min), num2str( (t2_max - t2_min) / 2), num2str(t2_max)}, ...
    {num2str(t2_min), num2str( (t2_max - t2_min) / 2), num2str(t2_max)}, ...
    {'0', '0.5', '1.0'}, ...
    };
%
%%% Figure design
%
par_t_mat  = [
    0  2  2  2  2  2  2
    1  0  2  2  2  2  2
    1  1  0  2  2  2  2
    1  1  1  0  2  2  2
    1  1  1  1  0  2  2
    1  1  1  1  1  0  2
    1  1  1  1  1  1  0];
%
par_x_mat  = [
    0  2  2  2  2  2  2
    2  0  3  3  3  3  3
    2  3  0  4  4  4  4
    2  3  4  0  5  5  5
    2  3  4  5  0  11 11
    2  3  4  5  11 0  12
    2  3  4  5  11 12 0];
%
par_y_mat  = [
    0  3  4  5  11 12 13
    3  0  4  5  11 12 13
    4  4  0  5  11 12 13
    5  5  5  0  11 12 13
    11 11 11 11 0  12 13
    12 12 12 12 12 0  13
    13 13 13 13 13 13 0];
%
[n_row, n_col] = size(par_t_mat);


%%% Plotting
%
figure_nr       = 788;
n_cols_figure   = 2;
journal         = 'mrm';
res_dpi         = 1000;
%
m_upper = 0;
ph      = 1;
phm     = 1;
m_lower = 1.2;
%
m_left  = 1.0;
pw      = 1;
pwm     = 1;
m_right = 0;


% Axis setup
fh = m_upper + n_row * ph + (n_row - 1) * phm + m_lower;
fw = m_left  + n_col * pw + (n_col - 1) * pwm + m_right;
%
m_upper = m_upper / fh;
ph      = ph / fh;
phm     = phm / fh;
m_lower = m_lower / fh;
%
m_left  = m_left / fw;
pw      = pw / fw;
pwm     = pwm / fw;
m_right = m_right / fw;
%
imsize     = my_imsize(n_cols_figure, fh / fw, journal);
%
fs_axis     = 6;
fs_label    = fs_axis;
lw          = 0.7;
%
fa = 0.4;
lw_scatter = 0.3;
ms      = 5;
ms_gt   = 5;
%
par_sc = [1 1 1e9 1e9 1 1 1 1 1 1 1e3 1e3 1];

%%%------------------------------------------------------------------
% PERFORM
%%%------------------------------------------------------------------

%%% Extract
%
my_struct   = load(fullfile(output_dir, 'smr_fig_si45_sim.mat'));
sim_mp      = my_struct.sim_mp;
n_plot      = numel(c_case_plot);



%%% Check for bias
%
c_mp_check = [1:5 11:12];
mp_fit = squeeze(mean(sim_mp(:,:,c_mp_check),2))';
tmp    = smr_mp(mp_list{1});
mp_gt  = repmat(tmp(c_mp_check)', [1 8]);
mp_bias = 100 * (mp_fit - mp_gt) ./ mp_gt
%        
mp_gt = [smr_mp(mp_list{1}) gt_p2] .* par_sc;  



%%% Plot
%
figure(figure_nr)
clf
set(gcf, 'color', 'w')

for c_row = 1:n_row
    %
    b_pos = 1 - m_upper - c_row * ph - (c_row - 1) * phm;
    
    for c_col = 1:n_col
        %
        l_pos = m_left + (c_col - 1) * (pw + pwm);
        
        
        if (par_t_mat(c_row, c_col) ~= 1), continue; end
        
        
        % Axes
        ax_l = l_pos;
        ax_b = b_pos;
        ax_w = pw;
        ax_h = ph;
        axes('position', [ax_l ax_b ax_w ax_h]);
        
        % Plot
        c_par_x = par_x_mat(c_row, c_col);
        c_par_y = par_y_mat(c_row, c_col);
        %        
        %
        for c_p = 1:n_plot
            c_case = c_case_plot(c_p);
            %
            vx = sim_mp(c_case, :, c_par_x) * par_sc(c_par_x);
            %
            if (c_par_y == 13)
                vy = smr_px2p2(squeeze(sim_mp(c_case,:,6:10)));
            else
                vy = sim_mp(c_case, :, c_par_y) * par_sc(c_par_y);
            end
            %
            scatter(vx, vy, ms, ...
                'markeredgecolor', ep_col_list{c_case}, 'markerfacecolor', ep_col_list{c_case}, ...
                'MarkerFaceAlpha', fa, 'linewidth', lw_scatter);
            hold on
        end
        
        % Plot ground truth
        plot(mp_gt(c_par_x), mp_gt(c_par_y), 'xk', 'markersize', ms_gt)
        %
        axis square;
        %
        set(gca, 'tickdir', 'out', 'fontsize', fs_axis, 'box', 'off', ...
            'linewidth', lw, 'fontweight', 'bold')
        %
        xl = par_range_list{c_par_x};
        yl = par_range_list{c_par_y};
        xlim(xl);
        ylim(yl);
        %
        set(gca, 'xtick',  linspace(xl(1), xl(2), numel(ticklbl_list{c_par_x})), 'xticklabel', ticklbl_list{c_par_x});
        set(gca, 'ytick',  linspace(yl(1), yl(2), numel(ticklbl_list{c_par_y})), 'yticklabel', ticklbl_list{c_par_y});
        
        xlabel(par_name_list{c_par_x}, 'fontsize', fs_label)
        ylabel(par_name_list{c_par_y}, 'fontsize', fs_label)
        
    end
end

%%% PRINT FIGURE TO FILE
output_filename = fullfile(output_dir, output_name);
save_current_fig_to_file(output_name, output_dir, imsize, res_dpi);
disp(['Wrote: ' output_filename]);
system(['open ' output_filename]);

end