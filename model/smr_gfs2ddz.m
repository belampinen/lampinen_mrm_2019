function ddz = smr_gfs2ddz(g, fs)
% function s = smr_gfs2ddz(g, fs)

fsp     = fs .* (1 + (g.^2 ./ (1 - g.^2))) ./ (fs + (g.^2 ./ (1 - g.^2)));
ddz     = fsp ./ (3 - 2 * fsp);

end