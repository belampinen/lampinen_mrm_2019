function g_pf = smr_optimize_penalize_g(g, g_limit)
%
do_penalize         = g > g_limit;
%
sc                  = 500000;
%
g_pf                = ((2 + sc * (g - g_limit)) .* do_penalize) + 1;
end