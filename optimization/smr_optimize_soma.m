function result = smr_optimize_soma(opt)

%-------------------------------------------------------------------------
% PREPARE
%-------------------------------------------------------------------------
%
%

%%% SOMA
soma_opt.fh_metric          = @(guess) smr_optimize_guess2metric(guess, opt);
soma_opt.fh_log             = @(msg)   smr_optimize_log(opt.fid, msg, ~opt.silence, 1, 0, 0);
soma_opt.fh_printguess      = @(guess) smr_optimize_printguess(guess, opt);
soma_opt.bounds             = opt.bounds;

% Changes to default SOMA settings
soma_opt.n_mig              = 200;
soma_opt.p                  = 0.2;


%-------------------------------------------------------------------------
% PERFORM
%-------------------------------------------------------------------------

v_pop        = cell(opt.n_run);
d_pop        = v_pop;
d_history    = v_pop;
%
for c_run = 1:opt.n_run
    %
    feval(soma_opt.fh_log, ['Run nr: ' num2str(c_run) ' -------------------']);
    feval(soma_opt.fh_log, ' ');
    %    
     [v_pop{c_run}, d_pop{c_run}, d_history{c_run}]  = soma_bl(soma_opt);
    %
    feval(soma_opt.fh_log, ' ');
    feval(soma_opt.fh_log, ' ');
end
%
result = {v_pop, d_pop, d_history, soma_opt};

end