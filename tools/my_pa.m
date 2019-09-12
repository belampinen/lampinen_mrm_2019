function s_pa = my_pa(s, xps)
s_ind   = unique(xps.s_ind);
s_pa    = zeros(size(s_ind));
for c = 1:numel(s_ind)
    ind = xps.s_ind == s_ind(c);
    s_pa(c) = mean(s(ind));
end
end