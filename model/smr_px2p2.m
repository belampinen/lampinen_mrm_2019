function p2 = smr_px2p2(px)
%
% px = [p20 Re(p21) Im(p21) Re(p22) Im(p22)]
%
if (size(px, 2) ~= 5), px = px'; end
%
p20     = px(:,1);
p2p1    = 1 * px(:,2) + 1i * px(:,3); 
p2p2    = 1 * px(:,4) + 1i * px(:,5); 
p2m1    = -conj(p2p1);
p2m2    = conj(p2p2);
%
p2m     = [p2m2 p2m1 p20 p2p1 p2p2];
%
%
N       = sqrt(5 / (4*pi));
p2      = sqrt(sum( (abs(p2m)).^2, 2) ) / N;
end