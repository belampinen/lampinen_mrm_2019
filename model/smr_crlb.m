function crlb = smr_crlb(xps, mp, snr)

%--------------------------------------------------------
% PREPARE
%--------------------------------------------------------

%%% Extract model parameters
s0      = mp(1);  
f       = mp(2);  
di_z    = mp(4);
dd_z    = mp(5);
t2_s    = mp(11);
t2_z    = mp(12);

%%% Derive quantities
[~, export] = smr_fit2data(mp, xps, 1);
%
a_s  = export.a_s;
a_z  = export.a_z;
A_s  = export.A_s;
A_z  = export.A_z;
C2_s = export.C2_s;
C2_z = export.C2_z;
y20  = export.y20;
y21  = export.y21;
y22  = export.y22;
pys  = export.pys;
Ads  = export.Ads;
Adz  = export.Adz;
At2s = export.At2s;
At2z = export.At2z;

% 
%--------------------------------------------------------
% PERFORM
%--------------------------------------------------------

%%% Model
% s   = s0 * [ f * Ads .* At2s + (1 - f) * Adz .* At2z]
% Ad  = exp(a) .* (1/2 * C0 + pi * C2 .* pys)
%

%%% Derivatives
%
% --- S0
% D1  = s / s0;
D1  = f * Ads .* At2s + (1 - f) * Adz .* At2z;

% --- f
D2  = s0 * (Ads .* At2s - Adz .* At2z);

% --- dx
a_s_di_s = -xps.b .* (1 - xps.b_delta * 1);
a_z_di_z = -xps.b .* (1 - xps.b_delta * dd_z);
a_z_dd_z =  xps.b .* xps.b_delta * di_z;
%
A_s_di_s = 3 * xps.b .* xps.b_delta * 1;
A_z_di_z = 3 * xps.b .* xps.b_delta * dd_z;
A_z_dd_z = 3 * xps.b .* xps.b_delta * di_z;
%
% Obtain the CX derivatives numerically
step            = 1e-4;
[C0_sf, C2_sf]  = smr_c0c2(A_s+step);
[C0_sb, C2_sb]  = smr_c0c2(A_s-step);
C0_s_p = (C0_sf - C0_sb) / (2 * step);
C2_s_p = (C2_sf - C2_sb) / (2 * step);
%
[C0_zf, C2_zf]  = smr_c0c2(A_z+step);
[C0_zb, C2_zb]  = smr_c0c2(A_z-step);
C0_z_p = (C0_zf - C0_zb) / (2 * step);
C2_z_p = (C2_zf - C2_zb) / (2 * step);
%
Ads_di_s = a_s_di_s .* Ads + exp(a_s) .* A_s_di_s .* (C0_s_p/2 + pi * C2_s_p .* pys);
Ads_di_z = a_z_di_z .* Adz + exp(a_z) .* A_z_di_z .* (C0_z_p/2 + pi * C2_z_p .* pys);
Ads_dd_z = a_z_dd_z .* Adz + exp(a_z) .* A_z_dd_z .* (C0_z_p/2 + pi * C2_z_p .* pys);
%
D3  = s0 * f          * Ads_di_s .* At2s;  % di_s
D4  = s0 * (1 - f)    * Ads_di_z .* At2z;  % di_z
D5  = s0 * (1 - f)    * Ads_dd_z .* At2z;  % dd_z

% --- p2x
c_f = s0 * pi * (f * exp(a_s) .* C2_s .* At2s + (1 - f) * exp(a_z) .* C2_z .* At2z);
%
D6  =           y20    .* c_f; % p20
D7  =  2 * real(y21)   .* c_f; % p21r
D8  = -2 * imag(y21)   .* c_f; % p21i
D9  =  2 * real(y22)   .* c_f; % p22r
D10 = -2 * imag(y22)   .* c_f; % p22i

% --- T2
D11 = s0 * f       * Ads .* (xps.te / t2_s^2) .* At2s;
D12 = s0 * (1 - f) * Adz .* (xps.te / t2_z^2) .* At2z;
%
%
%    s0     fs      di_s    di_z    dd_z    p20     p21r    p21i    p22r    p22i    t2_s    t2_z 
D = [D1     D2      D3      D4      D5      D6      D7      D8      D9      D10     D11     D12];


%%% Scale derivatives before inversion for better performance/conditioning
%
sc_vec = 1 ./ sqrt(sum(D.^2));
D = D.* repmat(sc_vec, [xps.n 1]);

%%% F-matrix elements from D-products summed across observations 
%
F = zeros(12, 12);
for c_row = 1:12
    for c_col = 1:12
        F(c_row, c_col) = sum(D(:,c_row) .* D(:,c_col));
    end
end
F = F * snr^2;  % s0 = 1 per default --> sigma_noise = 1 / snr

%%% CRLB
%
if (rcond(F) > eps)
    covar_matrix = inv(F);
    crlb = diag(covar_matrix);
else
    crlb = (zeros(size(mp)) + 1e9)'; % Something high but not inf
end
% Re-scale
crlb = crlb .* (sc_vec(:).^2);
end


