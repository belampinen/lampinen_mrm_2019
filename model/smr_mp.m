function mp = smr_mp(mp_name)

unit_to_SI = [1 1 1e-9 1e-9 1 1 1 1 1 1 1e-3 1e-3];

switch (mp_name)
    
    case 'a'  % represents (fairly coherent) WM
        % p2m = smr_od2p2m([0.1 0 1] / norm([0.1 0 1] ), 0.2) <-- od     
        p2m = [0.5576    0.0685    0.0001    0.0023    -0.0000];        
        %
        %       s0      fs      di_s    di_z    dd_z    p2m     t2_s    t2_z
        mp  =   [1      0.45    0.6     1.3     0.57    p2m     80      60];
        
    case 'b'  % represents DGM 
        % p2m = smr_od2p2m([0.1 0 1] / norm([0.1 0 1] ), 0.45) <-- od
        p2m =   [0.2110    0.0260    0.0000    0.0011   -0.0000];
        %       s0      fs      di_s    di_z    dd_z    p2m     t2_s    t2_z
        mp  =   [1      0.15    0.3     0.9     0.4     p2m     75      55];
        
    case 'c' % represents WML
        % p2m = smr_od2p2m([0.1 0 1] / norm([0.1 0 1] ), 0.35) <-- od
        p2m =   [0.2990    0.0368    0.0000    0.0015   -0.0000];
        %       s0      fs      di_s    di_z    dd_z    p2m     t2_s    t2_z
        mp  =   [1      0.4     0.6     1.7     0.4     p2m     80      150];
           
end

mp = mp .* unit_to_SI;

end