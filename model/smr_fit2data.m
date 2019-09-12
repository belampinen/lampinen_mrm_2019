function [s, export] = smr_fit2data(m, xps, do_export)
% function s = smr_fit2data(m, xps)
%
% s = s0 * [ f * Ads * At2s + (1 - f) * Adz * At2z]
%
if (nargin < 3), do_export = 0; end
export = [];

%%% Extract model parameters
s0      = m(1);
%
% Kernel
f       = m(2);
di_s    = m(3);
dd_s    = 1;
di_z    = m(4);
dd_z    = m(5);
%
% ODF
p20     = m(6);
p21r    = m(7);
p21i    = m(8);
p22r    = m(9);
p22i    = m(10);
%
% Relaxation
t2_s    = m(11);
t2_z    = m(12);
%
% Convert gradient vectors from Cartesian to Spherical
x       = xps.u(:,1);
y       = xps.u(:,2);
z       = xps.u(:,3);
[phi, theta] = cart2sph(x, y, z);
phi     = phi + pi;       % [-pi pi]  --> [0 2pi]
theta   = pi / 2 - theta; % elevation --> inclination

%%% Diffusion attenuation
%
a_s = -xps.b * di_s .* (1 - xps.b_delta * dd_s);
a_z = -xps.b * di_z .* (1 - xps.b_delta * dd_z);
%
A_s = 3 * xps.b .* xps.b_delta * di_s * dd_s;
A_z = 3 * xps.b .* xps.b_delta * di_z * dd_z;
%
[C0_s, C2_s] = smr_c0c2(A_s);
[C0_z, C2_z] = smr_c0c2(A_z);
%
y20 = smr_spha(2,  0, theta, phi);
y21 = smr_spha(2,  1, theta, phi);
y22 = smr_spha(2,  2, theta, phi);
pys    = ...
    p20 * y20 + ...
    p21r * ( 2) * real(y21) + ...
    p21i * (-2) * imag(y21) + ...
    p22r * ( 2) * real(y22) + ...
    p22i * (-2) * imag(y22); 
%
% first term is the pa (as in Eriksson 2015) - carries the signal weight
% second term is the deviation from the pa due to voxel-level anisotropy
% the px values relative to each other determines the shape of modulation
% the weight on all px changes the degree of directional modulation
Ads     = exp(a_s) .* (C0_s/2 + pi * C2_s .* pys);
Adz     = exp(a_z) .* (C0_z/2 + pi * C2_z .* pys); 

% Already accounted for:
% y00 = 1 / (2 * sqrt(pi))
% p00 = 1 / sqrt(pi);          
% p00 * y00 = 1 / (2 * pi)


%%% Relaxation attenuation
%
At2s = exp(-xps.te/t2_s);
At2z = exp(-xps.te/t2_z);


%%% Signal
s = s0 * (f * Ads .* At2s + (1 - f) * Adz .* At2z);

%%% Export
if (do_export)
    export.a_s  = a_s;
    export.a_z  = a_z;
    export.A_s  = A_s;
    export.A_z  = A_z;
    export.C0_s = C0_s;
    export.C0_z = C0_z;
    export.C2_s = C2_s;
    export.C2_z = C2_z;
    export.y20  = y20;
    export.y21  = y21;
    export.y22  = y22;
    export.pys  = pys;
    export.Ads  = Ads;
    export.Adz  = Adz;
    export.At2s = At2s;
    export.At2z = At2z;
end
end