function smr_fig_4(output_dir)

if (nargin < 1), output_dir = fullfile(pwd, 'figures', 'output'); end

if (0)
    smr_fig_4_sim;
end


%%%------------------------------------------------------------------
% PREPARE
%%%------------------------------------------------------------------

%
output_name = 'smr_fig_4.tiff';
%
mp_list = {'zero_od', 'mid_od', 'full_od', 'full_od_plus'}; 
n_mp        = numel(mp_list);
xticklbl    = {'Zero OD', 'Mid OD', 'Full OD', 'Full OD+'};
%
c_par_list   = [2 3 4 5 11 12 13]; % fs, dis, diz, ddz, T2s, T2z, p2
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
    ['\itD\rm\bf_{I;S} [' char(181) 'm^2/ms]'], ...
    ['\itD\rm\bf_{I;Z} [' char(181) 'm^2/ms]'], ...
    '\itD\rm\bf_{\Delta;Z}', ...
    [], ...
    [], ...
    [], ...
    [], ...
    [], ...
    'T_{2;S} [ms]', ...
    'T_{2;Z} [ms]', ...
    '\itp\rm\bf_2', ...
    'MSR'};
%
par_list_simple = {'S0', 'fs', 'dis', 'diz', 'ddz', [], [], [], [], [], 't2s', 't2z', 'p2', 'msr'};
%
par_sc    = [1 1 1e9 1e9 1 1 1 1 1 1 1e3 1e3 1 1];
par_do_gt = [0 1 1   1   1 0 0 0 0 0 1   1   1 0];   
%
msr_max = 300;
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
    [0 1], ...
    [0 msr_max], ... % msr
    };
par_tick_lbl = {...
    [], ...
    {'0.0', '0.2', '0.3', '0.4', '0.6', '0.8', '1.0'}, ...
    {'0.00', '0.25', '0.50', '0.75', '1.00', '1.25', '1.50'}, ...
    {'0.00', '0.25', '0.50', '0.75', '1.00', '1.25', '1.50'}, ...
    {'-0.50', '-0.25', '0.00', '0.25', '0.50', '0.75', '1.00'}, ...
    [], ...
    [], ...
    [], ...
    [], ...
    [], ...
    {'50', '60', '70', '80', '90', '100'}, ...
    {'50', '60', '70', '80', '90', '100'}, ...
    {'0.0', '0.2', '0.3', '0.4', '0.6', '0.8', '1.0'}, ...
    [], ...
    };

%%% Figure
%
n_cols_figure   = 2;
journal         = 'mrm';
res_dpi         = 600;
%
%
%
m_upper = 0.2;
ph      = 1;         % plot
phm     = 0.5;
m_lower = 0.5;
%
m_left  = 0.4;
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
p2_gt           = [0.8882 0.1442 0 0 0.45];
%
m_bs.ms     = 2;
m_bs.mec    = 'k';
x           = [1 2 3 4];





%%%------------------------------------------------------------------
% PERFORM
%%%------------------------------------------------------------------

%%% Extract boot
%
my_struct = load(fullfile(output_dir, 'smr_fig_4_sim.mat'));
%
pars      = zeros(n_mp, n_par, 200);
pars_sd   = zeros(n_mp, n_par);
pars_m    = pars_sd;
pars_gt   = pars_sd;
%
for c_mp = 1:n_mp
     %
     mpw     = smr_mp_watson(mp_list{c_mp});   
     %
     s       = squeeze(my_struct.sim_s(c_mp,:,:));
     s_fit   = squeeze(my_struct.sim_s_fit(c_mp,:,:));
     
    for c_p = 1:n_par                                
        %
        c_par   = c_par_list(c_p);
        c_par_w = c_par_list_w(c_p);
        
        %%% Ground truth
        % 
        if (par_do_gt(c_par))
            if (c_par == 13)
                pars_gt(c_mp, c_p) = p2_gt(c_mp);
            else
                pars_gt(c_mp, c_p) = mpw(c_par_w) * par_sc(c_par);
            end
        end
        
        %%% Statistics                       
        %
        if (c_par < 13)
            vals = squeeze(my_struct.sim_mp(c_mp, :, c_par)) * par_sc(c_par);            
        elseif (c_par == 13) % p2
            vals = smr_px2p2(squeeze(my_struct.sim_mp(c_mp, :, 6:10)));                
        elseif (c_par == 14) % MSR
            vals = mean((s - s_fit).^2, 2);
        end
        pars   (c_mp, c_p, :) = vals;
        pars_m (c_mp, c_p)    = mean(vals);
        pars_sd(c_mp, c_p)    = std(vals);                
    end
end



%%% Create figure
%
figure(675)
clf
set(gcf, 'color', 'w')
%
c_p         = 1;
%
for c_row = 1:n_row
    for c_col = 1:n_col
        
        if (c_p > n_par), continue; end
        
        c_par   = c_par_list(c_p);
        
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
        
        
        % Plot ground truth
        if (par_do_gt(c_par))
            for c_mp = 1:n_mp
                plot(linspace(c_mp-0.5, c_mp+0.5, 10), repmat(pars_gt(c_mp,c_p), [10 1]), '-.', ...
                    'color', col_gt, 'linewidth', lw_gt)
                hold on
            end
        end
        1;
        
        % Plot beeswarms
        for c_mp = 1:n_mp
            tmp = squeeze(pars(c_mp,c_p,:));
            my_plot_beeswarm(x(c_mp), {tmp}, [], [], m_bs, 0);
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
        disp([...
            par_list_simple{c_par} ' = ' num2str(pars_gt(1,c_p)) ...
            ', mean = ' num2str(pars_m(1, c_p)) ...
            ', sd = ' num2str(pars_sd(1,c_p))]);               
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