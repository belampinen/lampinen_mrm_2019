function str = smr_optimize_printguess(guess, opt)
%
ep  = smr_optimize_guess2ep(guess, opt);

%
n       = sprintf('\n');
n_shell = size(ep,1);
%
% Sort
[~, index] = sort(ep(:,1), 'ascend');
%
ew = 10;
%
str = [char(repmat(45, [1 (ew * 3)+4])) n];
str = [str my_fill('b', ew) my_fill('b_delta', ew) my_fill('te', ew) my_fill('n_dir', ew) n];
for c_shell = 1:n_shell
    for c_ep = 1:4
        str = [str my_fill(num2str(ep(index(c_shell), c_ep)), ew)];
    end
    str = [str n];
end
%
% Include some other info
[~, export] = smr_optimize_ep2metric(ep, opt);
%
str         = [str n 'Vw      = ' num2str(export.vw)];
str         = [str n 't_acq   = ' num2str(export.t_acq / 60, 2) ' min'];
str         = [str n 'g_max   = ' num2str(export.g_max*1e3, 2) ' mT/m, (g_pf = ' num2str(export.g_pf) ')'];
str         = [str n 'snr_min = ' num2str(export.snr_min, 2) ', (snr_pf = ' num2str(export.snr_pf) ')'];

end
