function mp = smr_mp_watson(mp_name)

unit_to_SI = [1 1 1e-9 1e-9 1 1 1e-3 1e-3];

switch (mp_name)

    case 'zero_od'
        
        % p2m = smr_od2p2m([0.1 0 1] / norm([0.1 0 1] ), 1)
        %    s0    f        dis     diz     ddz     od      t2s    t2z
        mp = [1    0.45     0.6     1.3     0.57    0.05    80     60];
        
    case 'mid_od'
        
        %    s0    f        dis     diz     ddz     od      t2s    t2z
        mp = [1    0.45     0.6     1.3     0.57    0.5     80     60];
        
    case 'full_od'
        
        %    s0    f        dis     diz     ddz     od      t2s    t2z
        mp = [1    0.45     0.6     1.3     0.57    1       80     60];
        
    case 'full_od_plus'
        
        %    s0    f        dis     diz     ddz     od      t2s    t2z
        mp = [1    0.45     0.95    0.95    0.57    1       70     70];
            
end

mp = mp .* unit_to_SI;

end