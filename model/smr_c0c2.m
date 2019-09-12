function [c0, c2] = smr_c0c2(A)
%
A_thresh    = 1e-6;
is_low      = abs(A) < A_thresh;
%
f_low       = @(A, n) (gamma(n + 1/2) / gamma(2*n + 3/2)) * (-A).^(n);
f_c0_high   = @(A)    sqrt(pi) * real(gammainc(A,1/2)./sqrt(A));
f_c0        = @(A)    is_low .* f_low(A,0) + ~is_low .* f_c0_high(A+eps);
%
f_c2_high   = @(A)    (3/2) ./ A .* ((1/2) * f_c0(A) - exp(-A)) - (1/2) * f_c0(A);
f_c2        = @(A)    is_low .* f_low(A,1) + ~is_low .* f_c2_high(A+eps);
%
c0          = f_c0(A);
c2          = f_c2(A);
end