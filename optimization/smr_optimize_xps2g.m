function g = smr_optimize_xps2g(xps, opt)
%
% Calculates the gradient amps required to produce a list of [b, b_delta, te]
%
%

gamma   = 2 * pi*42.6e6;

% b_delta penalty
% Empiricaly estimated at a Siemens Prisma 80 mT/m, slew 200 mT/m/s
f_k     = @(bd) 0.0500 + bd * 0.0273 + bd.^2 * 0.0491 + bd.^3 * (-0.2032 ) + bd.^4 *  0.3492;
k       = f_k(xps.b_delta);

% Assymetry penalty factor (increases with lower TE)
y1 = 0.4;   % TE < 60 ms 
y2 = 1;     % TE > 85 ms (no penalty)
x1 = 0.060; 
x2 = 0.085;
%
a  = (y2 - y1) / (x2 - x1);
b  = y1 - a * x1;
%
ind1 = xps.te <  x1;
ind2 = xps.te >= x2;
ind3 = ~ind1 & ~ind2;
%
k = k .* (ind1 * y1 + ind3 .* (a * xps.te + b) + ind2 * y2); 

% Calculate effective TE encoding, using unbalanced gradients
te_encoding     = max(xps.te - opt.te_nonenc, 0);

% Model from (Sj?lund 2015; Eq. 13)
% b = bd_factor * (gamma^2 * g^2 * te_encoding.^3) / 4;
g = (2 / gamma) * sqrt( xps.b ./ (k .* te_encoding.^3));


end





