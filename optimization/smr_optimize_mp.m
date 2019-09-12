function mp = smr_optimize_mp(mp_name)

%%% Choose intervals
switch (lower(mp_name))
    
    case {'a', 'b', 'c'}
        
        mp = smr_mp(mp_name);   
        
    case 'composite'
        
        mp_a  = smr_mp('a');  
        mp_b  = smr_mp('b');  
        mp_c  = smr_mp('c');  
        mp = [mp_a; mp_b; mp_c];
        
end

end