function [h1, h2] = my_plot_beeswarm(X, Y, fmt, n_step, m, do_lines, ol_index)

if (nargin < 3), fmt = []; end
if (nargin < 4), n_step = []; end
if (nargin < 5), m = []; end
if (nargin < 6), do_lines = []; end
if (nargin < 7), ol_index = {}; end

if (isempty(fmt)), fmt = '.'; end
if (isempty(do_lines)), do_lines = 1; end

if (iscell(Y))
    for c = 1:numel(Y)
        
        % Handle format
        if (iscell(fmt))
            fmt_in = fmt{c};
        else
            fmt_in = fmt;
        end
        
        % Handle outliers
        if (~isempty(ol_index))
            if (~iscell(ol_index)), ol_index = {ol_index}; end
            ol_tmp = ol_index{c};
        else
            ol_tmp = zeros(1, numel(Y{c}));
        end
        
        % Plot
        m_tmp = m;
        if (~isempty(m) && isfield(m, 'mfc') && iscell(m.mfc))
            m_tmp.mfc = m.mfc{c};
        end
        [h1, h2] = my_plot_beeswarm(X(c), Y{c}, fmt_in, n_step, m_tmp, do_lines, ol_tmp); hold on;
    end
    return;
end

if (isempty(n_step))
    n_step = max(round(numel(Y) / 5), 1);
end

if (isempty(Y))
    h1 = [];
    h2 = [];
    return;
end

dX = zeros(size(Y));
dx_multiplier = 1;
if (~isempty(m) && isfield(m, 'dx_multiplier')), dx_multiplier = m.dx_multiplier; end

% Make beeswarm plot
mY = [min(Y) max(Y)];
dY = diff(mY) / n_step;

for y = (mY(1) - dY):dY:(mY(2) + dY)
    
    ind = find( (Y > y) & (Y <= (y + dY)));
    
    if (numel(ind) <= 1), continue; end
    
    f_w = @(n)(  (n - 1) * log(n) );
    w = f_w(numel(ind));
    
    [~,ind_sort] = sort(Y(ind));
    ind_sort = (ind_sort - 1) .* ((mod(ind_sort,2) * 2) - 1);
    ind_sort = ind_sort / max(abs(ind_sort)) / 2;
    
    dX(ind) = ind_sort * w;
    
end

if (numel(dX) ~= 1)
    dX = dX / std(dX) / 10;
end

if (isnan(dX(1)))
    %dX = [0 0.01];
    dX = zeros(size(dX));
end
dX = dX * dx_multiplier;

h1 = plot(X + dX((1 - ol_index) > 0), Y((1 - ol_index) > 0), fmt);
hold on
h2 = plot(X + dX(ol_index > 0), Y(ol_index > 0), fmt);
hold on;

if (do_lines)
    plot(X + [-1 1] * max(abs(dX)) * 1.6, mean(Y) + [0 0], 'k-', 'LineWidth', 2);
    plot(X + [-1 1] * max(abs(dX)), mean(Y) - std(Y) + [0 0], 'k-', 'LineWidth', 2);
    plot(X + [-1 1] * max(abs(dX)), mean(Y) + std(Y) + [0 0], 'k-', 'LineWidth', 2);
    plot(X + [0 0], mean(Y) + [-1 1] * std(Y), 'k-', 'LineWidth', 2);
end
if (~isempty(m) && isfield(m, 'mfc') && ~isempty(m.mfc))
    set(h1, 'MarkerFaceColor', m.mfc);
    set(h2, 'MarkerFaceColor', 'red');
end
if (~isempty(m) && isfield(m, 'mec') && ~isempty(m.mec))
    set(h1, 'MarkerEdgeColor', m.mec);
end
if (~isempty(m) && isfield(m, 'ms')  && ~isempty(m.ms))
    set(h1, 'MarkerSize', m.ms);
    set(h2, 'MarkerSize', m.ms);
end
end

