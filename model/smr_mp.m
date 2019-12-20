function mp = smr_mp(mp_name)

unit_to_SI = [1 1 1e-9 1e-9 1 1 1 1 1 1 1e-3 1e-3];

switch (mp_name)
    
    case 'a'  % represents (fairly coherent) WM
        %
        p2m = [0.2788    0.0343    0.0000    0.0012         0];
        % smr_plot_p2m(p2m)
        %
        %       s0      fs      di_s    di_z    dd_z    p2m     t2_s    t2_z
        mp  =   [1      0.45    0.6     1.3     0.57    p2m     80      60];
        
    case 'b'  % represents DGM 
        %
        p2m = [0.1055    0.0130         0    0.0005         0];       
        % smr_plot_p2m(p2m)
        %
        %       s0      fs      di_s    di_z    dd_z    p2m     t2_s    t2_z
        mp  =   [1      0.15    0.3     0.9     0.4     p2m     75      55];
        
    case 'c' % represents WML
        %
        p2m = [0.1495    0.0184         0    0.0007         0];
        % smr_plot_p2m(p2m)
        %
        %       s0      fs      di_s    di_z    dd_z    p2m     t2_s    t2_z
        mp  =   [1      0.4     0.6     1.7     0.4     p2m     80      150];
                           
    case 'a_full_od'
        %
        p2m = [0         0         0         0              0];
        % smr_plot_p2m(p2m)
        %
        %       s0      fs      di_s    di_z    dd_z    p2m     t2_s    t2_z
        mp  =   [1      0.45    0.6     1.3     0.57    p2m     80      60];
        
    case 'a_mid_od'
        %
        p2m = [0.0896    0.0111         0    0.0008   -0.0000];
        % smr_plot_p2m(p2m)
        %
        %       s0      fs      di_s    di_z    dd_z    p2m     t2_s    t2_z
        mp  =   [1      0.45    0.6     1.3     0.57    p2m     80      60];
        
        
end

mp = mp .* unit_to_SI;

end