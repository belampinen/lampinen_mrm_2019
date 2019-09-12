function p2m = smr_od2p2m(my, od, do_plot)

if (nargin < 3), do_plot = 1; end

my = my(:);

%%% Define experiment and model
xps             = my_ep2xps(smr_ep('in_vivo'));
%
opt             = smr_opt;
opt.smr.n_rep   = 5;
udirs           = uvec_elstat_500;
%
mp_watson       = smr_mp_watson('mid_od');
mp_watson(6)    = od;


%%% Generate signal with my and od
s               = smr_fit2data_watson(mp_watson, xps, my, udirs);

%%% Fit signal using the smr model
mp_smr          = smr_data2fit(s, xps, opt);
p2m             = mp_smr(6:10);
%
s_fit           = smr_fit2data(mp_smr, xps);


%%% Plot
if (do_plot)
    figure(627)
    clf
    set(gcf, 'color', 'w')
    %
    subplot(1,3,1)
    smr_plot_watson(od, my, 0)
    grid on
    title(['OD = ' num2str(od)])
    %
    subplot(1,3,2)
    smr_plot_p2m(p2m, 0)
    grid on
    title(['p2 = ' num2str(smr_p2m2p2(p2m))]);
%
    %
    subplot(1,3,3)
    plot(s, '.k');
    hold on
    plot(s_fit, 'or');
    title(['SS = ' num2str(sum((s - s_fit).^2))])
    set(gca, 'box', 'off');
    axis square
    %
    set(gcf, 'pos', [21   520   960   281])    
end
end