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
i2_s = export.i2_s;
i2_z = export.i2_z;
y20  = export.y20;
y21  = export.y21;
y22  = export.y22;
p2m_y2m_sum  = export.p2m_y2m_sum;
Ads  = export.Ads;
Adz  = export.Adz;
At2s = export.At2s;
At2z = export.At2z;
%
p00  = 1 / sqrt(4 * pi);
y00  = 1 / sqrt(4 * pi);

% 
%--------------------------------------------------------
% PERFORM
%--------------------------------------------------------

%%% Model
% s   = s0 * [ f * Ads .* At2s + (1 - f) * Adz .* At2z]
% Adx = exp(a_x) .* (i0_x + 4pi * i2_x .* p2m_y2m_sum);
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
% Obtain the il_dx derivatives numerically
step            = 1e-4;
[i0_sf, i2_sf]  = smr_i0i2(A_s+step);
[i0_sb, i2_sb]  = smr_i0i2(A_s-step);
i0_s_p = (i0_sf - i0_sb) / (2 * step);
i2_s_p = (i2_sf - i2_sb) / (2 * step);
%
[i0_zf, i2_zf]  = smr_i0i2(A_z+step);
[i0_zb, i2_zb]  = smr_i0i2(A_z-step);
i0_z_p = (i0_zf - i0_zb) / (2 * step);
i2_z_p = (i2_zf - i2_zb) / (2 * step);
%
Ads_di_s = a_s_di_s .* Ads + exp(a_s) .* A_s_di_s .* (i0_s_p * p00 * y00 * sqrt(4*pi) + i2_s_p .* p2m_y2m_sum * sqrt(4*pi/5));
Adz_di_z = a_z_di_z .* Adz + exp(a_z) .* A_z_di_z .* (i0_z_p * p00 * y00 * sqrt(4*pi) + i2_z_p .* p2m_y2m_sum * sqrt(4*pi/5));
Adz_dd_z = a_z_dd_z .* Adz + exp(a_z) .* A_z_dd_z .* (i0_z_p * p00 * y00 * sqrt(4*pi) + i2_z_p .* p2m_y2m_sum * sqrt(4*pi/5));
%
D3  = s0 * f          * Ads_di_s .* At2s;  % di_s
D4  = s0 * (1 - f)    * Adz_di_z .* At2z;  % di_z
D5  = s0 * (1 - f)    * Adz_dd_z .* At2z;  % dd_z

% --- p2x
c_f = s0 * (f * exp(a_s) .* i2_s * sqrt(4*pi/5) .* At2s + (1 - f) * exp(a_z) .* i2_z * sqrt(4*pi/5) .* At2z);

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


