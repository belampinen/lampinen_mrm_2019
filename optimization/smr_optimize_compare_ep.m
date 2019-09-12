function smr_compare_ep(ep1, ep2, opt)

if (nargin < 3), opt = smr_optimize_opt_derive(smr_optimize_opt); end

%
[of1, export1] = smr_optimize_ep2metric(ep1, opt);
[of2, export2] = smr_optimize_ep2metric(ep2, opt);
1;

%%% Print comparison
%
n   = sprintf('\n');
%
str         = [n '----------------------------------'];
str         = [str n 'tacq   : ' num2str(export1.t_acq / 60, 2) '  / ' num2str(export2.t_acq / 60, 2) '  [min]'];
str         = [str n 'Vw     : ' num2str(export1.vw, 2) '  / ' num2str(export2.vw, 2)];
str         = [str n 'snr_min: ' num2str(export1.snr_min, 2)  ' / ' num2str(export2.snr_min, 2)];
str         = [str n 'snr_pf : ' num2str(export1.snr_pf, 2)  ' / ' num2str(export2.snr_pf, 2)];
str         = [str n 'gmax   : ' num2str(export1.g_max, 2)  ' / ' num2str(export2.g_max, 2)];
str         = [str n 'gmax_pf: ' num2str(export1.g_pf, 2)  ' / ' num2str(export2.g_pf, 2)];
str         = [str n 'OF     : ' num2str(of1, 2) '  / ' num2str(of2, 2)];
str         = [str n '----------------------------------' n];
%
fprintf(str);




end