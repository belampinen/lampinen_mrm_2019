function uvec = my_uvec_platonic(n_dir)

switch (n_dir)
    
    case 1
        uvec = [0 0 1];
        
    case 6
        uvec = uvec_icosa;
        
    case 10
        uvec = uvec_dodeca;
        
    case 15
        uvec = uvec_icosaedge;
        
    case 16
        uvec = [...
            uvec_icosa
            uvec_dodeca];
        
    case 21
        uvec = [...
            uvec_icosa
            uvec_icosaedge];
        
    case 25
        uvec = [...
            uvec_dodeca;
            uvec_icosaedge];
        
    case 30
        uvec = uvec_icosadodecaedge;
        
    case 40
        uvec = [...
            uvec_dodeca;
            uvec_icosadodecaedge];
        
    case 36
        uvec = [...
            uvec_icosa
            uvec_icosadodecaedge];
        
    case 45
        uvec = [...
            uvec_icosaedge;
            uvec_icosadodecaedge];
        
    otherwise
        error(['my_uvec_platonic not defined for ' num2str(n_dir) ' dirs']);
end


end