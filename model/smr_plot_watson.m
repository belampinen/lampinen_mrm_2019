function smr_plot_watson(od, my, do_create_figure)

if (nargin < 3), do_create_figure = 1; end
if (nargin < 2), my = [0 0 1]'; end
if (nargin < 1), od = [1 0.84 0.5 0.16 0.04]; end

if (do_create_figure)
    figure(3347)
    clf
    set(gcf, 'color', 'w')
end

% Obtain directions
n_dirs = 100;
r = 1;
theta   = linspace(0, 2*pi, n_dirs)';
phi     = linspace(-pi/2, pi/2, n_dirs)';

z_cos   = kron(r * cos(phi), ones(n_dirs, 1));
x       = repmat(r * cos(theta), [n_dirs 1]) .* z_cos;
y       = repmat(r * sin(theta), [n_dirs 1]) .* z_cos;
z       = kron(r * sin(phi), ones(n_dirs, 1));
udirs = [x y z];

% Plot distribution for different OD
for c_od = 1:numel(od)
    
    kappa = 1 ./ tan((pi / 2) * max(od(c_od), eps));
    kappa(kappa > 700) = 700;  %

    
    % Weights
    w       = exp(kappa * (udirs * my).^2);
    w       = w / max(w);
    wm      = repmat(w, [1 3]);
    wdirs   = udirs .* wm;    

    % Plot distribution
    plot3(wdirs(:,1), wdirs(:,2), wdirs(:,3), '.');
    my_plot = max(w(:)) * my;
    hold on
end

p = plot3(my_plot(1), my_plot(2), my_plot(3), 'r.', 'markersize', 20);
xlabel('X');
ylabel('Y');
zlabel('Z');
set(gcf, 'color', 'white');
axis equal
end

