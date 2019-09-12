function [v_pop, d_pop, d_history] = soma_bl(opt)
%
% Implementation of the all-to-one SOMA

%-------------------------------------------------------------------------
% PREPARE
%-------------------------------------------------------------------------

%%% Init
%
if(~isfield(opt, 'n_mig')),         opt.n_mig = 100; end
if(~isfield(opt, 'step_length')),   opt.step_length = 0.21; end   
if(~isfield(opt, 'path_length')),   opt.path_length = 2.1; end     
if(~isfield(opt, 'p')),             opt.p = 0.1; end                 % probability of element mutation
if(~isfield(opt, 'min_decrease')),  opt.min_decrease = 0.00001; end  % secondary termination condition
%
%
n_dim           = size(opt.bounds, 1);  % number of variables to optize
n_ind           = 2 * n_dim;            % number of individuals

%%% Print some information
%
feval(opt.fh_log, '-------------------------------');
feval(opt.fh_log, 'Running SOMA all-to-one with settings:');
feval(opt.fh_log, ['n_mig: ' num2str(opt.n_mig)]);
feval(opt.fh_log, ['n_ind: ' num2str(n_ind)]);
feval(opt.fh_log, ['step_length: ' num2str(opt.step_length)]);
feval(opt.fh_log, ['path_length: ' num2str(opt.path_length)]);
feval(opt.fh_log, ['p: ' num2str(opt.p)]);
feval(opt.fh_log, ['min_decrease: ' num2str(opt.min_decrease)]);
feval(opt.fh_log, '-------------------------------');
feval(opt.fh_log, ' ');


%-------------------------------------------------------------------------
% PERFORM
%-------------------------------------------------------------------------

%%% Initialize populations randomly from uniform distributions within bounds
%
while (1)
    v_pop   = zeros(n_dim, n_ind); % individuals (v)
    d_pop   = zeros(n_ind, 1);     % distances   (d)
    for c_ind = 1:n_ind        
        v_pop(:,c_ind)  = opt.bounds(:,1) + rand(size(opt.bounds(:,1))) .* (opt.bounds(:,2) - opt.bounds(:,1));
        d_pop(c_ind)    = feval(opt.fh_metric, v_pop(:,c_ind));
    end
    if (any(d_pop ~= Inf)), break; end
end


%%% Perform migrations
%
d_history    = zeros(opt.n_mig,1);
%
for c_mig = 1:opt.n_mig
    
    % Identify current leader (vl) = the individual with the lowest d
    c_ind_leader = find(d_pop == min(d_pop));   
    c_ind_leader = c_ind_leader(1);
    v_l          = v_pop(:,c_ind_leader);
    
    % Migrate non-leader individuals towards the leader by checking
    % by adjusting randomly chosen vector elements to become 
    % more similar to those of the leader
    for c_ind = 1:n_ind
        
        % Skip the leader
        if c_ind == c_ind_leader, continue; end
        
        % Store the individual's start position
        v_0 = v_pop(:,c_ind);
        
        % Explore a range of mutated versions
        for step_size = 0:opt.step_length:opt.path_length             
            % step_size = step number (m) x step length (sl)
                                    
            % 1) Generate vector designating which elements to take steps in,
            % each element has probability opt.p of stepping
            % Inlcuding step_size = 0 leaves option 
            x = zeros(size(v_0));
            while (sum(x) == 0)
                for c_dim = 1:n_dim
                    if rand < opt.p
                        x(c_dim) = 1;
                    end
                end
            end
                      
            % 2) Each element in mutation moves 'step_size' 
            % towards leader (with probability opt.p)
            v_m = v_0 + step_size * (v_l - v_0) .* x;
            
            % 3) Set values outside boundaries to random value
            %
            v_rand  = opt.bounds(:,1) + rand(size(opt.bounds(:,1))) .* (opt.bounds(:,2) - opt.bounds(:,1));
            ind     = (v_m < opt.bounds(:,1)) | (v_m > opt.bounds(:,2));
            v_m(ind) = v_rand(ind);
            
            % 4) Update position if better than previous
            d_m = feval(opt.fh_metric, v_m);
            if (d_m < d_pop(c_ind))
                d_pop(c_ind)   = d_m;   
                v_pop(:,c_ind) = v_m;
            end
        end
    end        
    
    
    % Identify current leader
    c_ind_leader = find(d_pop == min(d_pop)); 
    c_ind_leader = c_ind_leader(1);
  
    %%% LOGGING
    feval(opt.fh_log, ' ');
    feval(opt.fh_log, ['migration = '    num2str(c_mig)]);
    feval(opt.fh_log, 'Best guess:');
    feval(opt.fh_log, feval(opt.fh_printguess, v_pop(:,c_ind_leader)));
    feval(opt.fh_log, ['metric: '        num2str(d_pop(c_ind_leader))]);
        
    %%% Check secondary termination condition
    %
    d_history(c_mig) = d_pop(c_ind_leader);
    %
    n_mig_check = 5;
    if (c_mig > n_mig_check - 1)                
        decrease = abs(d_history(c_mig) - d_history(c_mig-n_mig_check+1));
        if (decrease < opt.min_decrease)            
            feval(opt.fh_log, 'Secondary termination condition reached');
            return
        end
    end
end
feval(opt.fh_log, 'Primary termination condition reached');
end