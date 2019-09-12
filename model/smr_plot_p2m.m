function p2 = smr_plot_p2m(px, do_create_figure)

if (nargin < 2), do_create_figure = 1; end

if (do_create_figure)
    figure(74)
    clf
    set(gcf, 'color', 'w')
end

% Obtain angles from directions
udirs = uvec_elstat_500;
%
[phi, theta] = cart2sph(udirs(:,1), udirs(:,2), udirs(:,3));
phi     = phi + pi;       % [-pi pi]  --> [0 2pi]
theta   = pi / 2 - theta; % elevation --> inclination

% Construct complex coefficients of ODF
p00     = 1 / sqrt(pi);
p20     = px(1);
p21r    = px(2);
p21i    = px(3);
p22r    = px(4);
p22i    = px(5);
%
y00 = smr_spha(0,  0, theta, phi);
y20 = smr_spha(2,  0, theta, phi);
y21 = smr_spha(2,  1, theta, phi);
y22 = smr_spha(2,  2, theta, phi);
odf    = ...
    p00 * y00 + ...
    p20 * y20 + ...
    p21r * ( 2) * real(y21) + ...
    p21i * (-2) * imag(y21) + ...
    p22r * ( 2) * real(y22) + ...
    p22i * (-2) * imag(y22); 
%
v = udirs .* repmat(odf, [1 3]);
plot3(v(:,1), v(:,2), v(:,3),  '.k');
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

max_val = max(abs([get(gca, 'xlim') get(gca, 'ylim') get(gca, 'zlim')]));
set(gca, 'xlim', [-max_val max_val], 'ylim', [-max_val max_val], 'zlim', [-max_val max_val])

set(gcf, 'pos', [ 21   435   441   366]);
%
p2 = smr_p2m2p2([p20 p21r p21i p22r p22i]);
%
title(['p2 = ' num2str(p2)])
end