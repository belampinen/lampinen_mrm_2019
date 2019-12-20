function g = smr_ddzfs2g(ddz, fs)
% function m = smr_ddzfs2g(ddz, fs)

t_lb = 0;
t_ub = 1;

    function ddz_fit = my_1d_fit2data(t,varargin)
        ddz_fit = smr_gfs2ddz(t, fs);
        1;
    end


lsqopts = optimoptions('lsqcurvefit', 'display', 'off','MaxFunEvals',1e4);

guess = rand(size(fs));
%
g = lsqcurvefit(@(t,varargin) my_1d_fit2data(t), guess, [], ...
    ddz, t_lb, t_ub, lsqopts);

end