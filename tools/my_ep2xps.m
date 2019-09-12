function xps = my_ep2xps(ep)

%%% Init the xps
%
xps.n       = sum(ep(:,4));
xps.b       = zeros(xps.n, 1);
xps.b_delta = xps.b;
xps.te      = xps.b;
xps.u       = zeros(xps.n, 3);
xps.s_ind   = xps.b;

%%% Fill out the xps
n_shell     = size(ep, 1);
%
c_start     = 1;
for c_shell = 1:n_shell

    n_dir = ep(c_shell, 4);
    c_stop = c_start + n_dir - 1;
    %
    xps.b      (c_start:c_stop)     = ep(c_shell, 1) * 1e9;
    xps.b_delta(c_start:c_stop)     = ep(c_shell, 2);
    xps.te     (c_start:c_stop)     = ep(c_shell, 3) * 1e-3;
    xps.u      (c_start:c_stop, :)  = uvec_platonic(n_dir);
    xps.s_ind  (c_start:c_stop)     = repmat(c_shell, [n_dir 1]); 
    %
    c_start = c_stop + 1;
end

end