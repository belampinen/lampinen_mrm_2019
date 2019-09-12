function vals_out = smr_optimize_discretize_ep(vals_in, ep_name, opt)

% Default, do nothing
vals_out = vals_in;

% Discetize, depending on ep_name
switch (ep_name)
    
    case 'b'
        
        space = opt.b_min:opt.b_step:opt.b_max;
        
    case 'b_delta'
        
        space = opt.b_delta_min:opt.b_delta_step:opt.b_delta_max;
        
    case 'te'
        
        space = opt.te_min:opt.te_step:opt.te_max;
        
    case 'ndir_scale'

        space = opt.ndir_scale_min:opt.ndir_scale_step:opt.ndir_scale_max;        
        
    otherwise
        space = [];
        
        
end
if (numel(space) > 1)
    vals_out = my_discretize(vals_in, space)';
end

end

%-------------------------------------
function vp = my_discretize(v, space)
% Maps vector v to values from vector 'space'

v       = v(:);
space   = space(:);

vgrid      = repmat(v', [numel(space) 1]);
spacegrid  = repmat(space, [1 numel(v)]);
[~, index]  = sort(abs(vgrid - spacegrid), 'ascend');
index = index(1,:);

vp = space(index);


end