function smr_optimize_master(output_dir, test_nr, n_worker, do_parallel)
%
if (nargin < 1), output_dir = fullfile(pwd, 'optimization', 'output'); end
if (nargin < 2), test_nr     = 7; end 
if (nargin < 3), n_worker    = 1; end
if (nargin < 4), do_parallel = 0; end
%
if (numel(test_nr) > 1)
    if (n_worker > 1)
        parpool(n_worker);
        parfor c_test = 1:numel(test_nr)
            smr_optimize_master(test_nr(c_test), 1, 1);
        end
        my_parpool(-1);
    else
        for c_test = 1:numel(test_nr)
            smr_optimize_master(test_nr(c_test));
        end
    end
    return;
end


%--------------------------------------------------
% PREPARE
%--------------------------------------------------

% Misc
%
do_rerun            = 1;


% Standard settings
%
opt                         = smr_optimize_opt;
opt.silence                 = 0;
opt.do_parallel             = do_parallel;
opt.discrete_ndir           = [1 6 10 15 30 45];
opt.n_shell                 = 10;

% Test-specific settings
%
switch test_nr
        
    %%% SHELL-SERIES
    case 1 % 4 shells
        test_name           = 'opt_4_shell';
        opt.n_shell         = 4;
        
    case 2 % 5 shells
        test_name           = 'opt_5_shell';
        opt.n_shell         = 5;
        
    case 3 % 6 shells
        test_name           = 'opt_6_shell';
        opt.n_shell         = 6;
        
    case 4 % 7 shells
        test_name           = 'opt_7_shell';
        opt.n_shell         = 7;
        
    case 5 % 8 shells
        test_name           = 'opt_8_shell';
        opt.n_shell         = 8;
        
    case 6 % 9 shells
        test_name           = 'opt_9_shell';
        opt.n_shell         = 9;
        
    case 7 % 10 shells / free b_delta
        test_name           = 'opt_10_shell';
        opt.n_shell         = 10;
        
        
        %%% B-delta series
    case 8 % lte
        test_name           = 'opt_lte';
        opt.b_delta_min     = 1;
        
    case 9 % dde
        test_name           = 'opt_dde';
        opt.discrete_b_delta = [-0.5 1];
        
    case 10 % ste
        test_name           = 'opt_ste';
        opt.discrete_b_delta = [0 1];
        
        
        %%% B_max series
        % g_max = 75
    case 11 % b2000
        test_name           = 'opt_b2000_g75';
        opt.b_max           = 2;
        
    case 12 % b2500
        test_name           = 'opt_b2500_g75';
        opt.b_max           = 2.5;
        
    case 13 % b3000
        test_name           = 'opt_b3000_g75';
        opt.b_max           = 3;
        
    case 14 % b4000
        test_name           = 'opt_b4000_g75';
        opt.b_max           = 4;
        
    case 15 % b5000
        test_name           = 'opt_b5000_g75';
        opt.b_max           = 5;
        
    case 16 % b6000
        test_name           = 'opt_b6000_g75';
        opt.b_max           = 6;
        
    case 17 % b7000
        test_name           = 'opt_b7000_g75';
        opt.b_max           = 7;
        
        % g_max = 200
    case 18 % b2000
        test_name           = 'opt_b2000_g200';
        opt.b_max           = 2.0;
        opt.g_max           = 0.2;
        
    case 19 % b2500
        test_name           = 'opt_b2500_g200';
        opt.b_max           = 2.5;
        opt.g_max           = 0.2;
        
    case 20 % b3000
        test_name           = 'opt_b3000_g200';
        opt.b_max           = 3;
        opt.g_max           = 0.2;
        
    case 21 % b4000
        test_name           = 'opt_b4000_g200';
        opt.b_max           = 4;
        opt.g_max           = 0.2;
        
    case 22 % b5000
        test_name           = 'opt_b5000_g200';
        opt.b_max           = 5;
        opt.g_max           = 0.2;
                
    case 23 % b6000
        test_name           = 'opt_b6000_g200';
        opt.b_max           = 6;
        opt.g_max           = 0.2;
                
    case 24 % b7000
        test_name           = 'opt_b7000_g200';
        opt.b_max           = 7;
        opt.g_max           = 0.2;
                                
        % g_max = 1e6 (~inf)  
    case 25 % b2000
        test_name           = 'opt_b2000_ginf';
        opt.b_max           = 2.0;
        opt.g_max           = 1e6;
        
    case 26 % b2500
        test_name           = 'opt_b2500_ginf';
        opt.b_max           = 2.5;
        opt.g_max           = 1e6;
        
    case 27 % b3000
        test_name           = 'opt_b3000_ginf';
        opt.b_max           = 3;
        opt.g_max           = 1e6;
        
    case 28 % b4000
        test_name           = 'opt_b4000_ginf';
        opt.b_max           = 4;
        opt.g_max           = 1e6;
        
    case 29 % b5000
        test_name           = 'opt_b5000_ginf';
        opt.b_max           = 5;
        opt.g_max           = 1e6;
                
    case 30 % b6000
        test_name           = 'opt_b6000_ginf';
        opt.b_max           = 6;
        opt.g_max           = 1e6;
                
    case 31 % b7000
        test_name           = 'opt_b7000_ginf';
        opt.b_max           = 7;
        opt.g_max           = 1e6;
        
                
        %%% TEMIN SERIES
        %  g_max = 75        
    case 32 % TEmin 10
        test_name           = 'opt_te10_g75';
        opt.te_min          = 10;
        
    case 33 % TEmin 20
        test_name           = 'opt_te20_g75';
        opt.te_min          = 20;
        
    case 34 % TEmin 30
        test_name           = 'opt_te30_g75';
        opt.te_min          = 30;
        
    case 35 % TEmin 40
        test_name           = 'opt_te40_g75';
        opt.te_min          = 40;
        
    case 36 % TEmin 50
        test_name           = 'opt_te50_g75';
        opt.te_min          = 50;
        
    case 37 % TEmin 60
        test_name           = 'opt_te60_g75';
        opt.te_min          = 60;
        
    case 38 % TEmin 70
        test_name           = 'opt_te70_g75';
        opt.te_min          = 70;
        
    case 39 % TEmin 80
        test_name           = 'opt_te80_g75';
        opt.te_min          = 80;
        
    case 40 % TEmin 90
        test_name           = 'opt_te90_g75';
        opt.te_min          = 90;
        
    case 41 % TEmin 100
        test_name           = 'opt_te100_g75';
        opt.te_min          = 100;
                
        % g_max = 200
    case 42 % TEmin 10
        test_name           = 'opt_te10_g200';
        opt.te_min          = 10;
        opt.g_max           = 0.2;
        
    case 43 % TEmin 20
        test_name           = 'opt_te20_g200';
        opt.te_min          = 20;
        opt.g_max           = 0.2;
        
    case 44 % TEmin 30
        test_name           = 'opt_te30_g200';
        opt.te_min          = 30;
        opt.g_max           = 0.2;
        
    case 45 % TEmin 40
        test_name           = 'opt_te40_g200';
        opt.te_min          = 40;
        opt.g_max           = 0.2;
        
    case 46 % TEmin 50
        test_name           = 'opt_te50_g200';
        opt.te_min          = 50;
        opt.g_max           = 0.2;
        
    case 47 % TEmin 60
        test_name           = 'opt_te60_g200';
        opt.te_min          = 60;
        opt.g_max           = 0.2;
        
    case 48 % TEmin 70
        test_name           = 'opt_te70_g200';
        opt.te_min          = 70;
        opt.g_max           = 0.2;
        
    case 49 % TEmin 80
        test_name           = 'opt_te80_g200';
        opt.te_min          = 80;
        opt.g_max           = 0.2;
        
    case 50 % TEmin 90
        test_name           = 'opt_te90_g200';
        opt.te_min          = 90;
        opt.g_max           = 0.2;
        
    case 51 % TEmin 100
        test_name           = 'opt_te100_g200';
        opt.te_min          = 100;
        opt.g_max           = 0.2;
                        
        % g_max = 1e6 (~inf)
    case 52 % TEmin 10
        test_name           = 'opt_te10_ginf';
        opt.te_min          = 10;
        opt.g_max           = 1e6;
        
    case 53 % TEmin 20
        test_name           = 'opt_te20_ginf';
        opt.te_min          = 20;
        opt.g_max           = 1e6;
        
    case 54 % TEmin 30
        test_name           = 'opt_te30_ginf';
        opt.te_min          = 30;
        opt.g_max           = 1e6;
        
    case 55 % TEmin 40
        test_name           = 'opt_te40_ginf';
        opt.te_min          = 40;
        opt.g_max           = 1e6;
        
    case 56 % TEmin 50
        test_name           = 'opt_te50_ginf';
        opt.te_min          = 50;
        opt.g_max           = 1e6;
        
    case 57 % TEmin 60
        test_name           = 'opt_te60_ginf';
        opt.te_min          = 60;
        opt.g_max           = 1e6;
        
    case 58 % TEmin 70
        test_name           = 'opt_te70_ginf';
        opt.te_min          = 70;
        opt.g_max           = 1e6;
        
    case 59 % TEmin 80
        test_name           = 'opt_te80_ginf';
        opt.te_min          = 80;
        opt.g_max           = 1e6;
        
    case 60 % TEmin 90
        test_name           = 'opt_te90_ginf';
        opt.te_min          = 90;
        opt.g_max           = 1e6;
        
    case 61 % TEmin 100
        test_name           = 'opt_te100_ginf';
        opt.te_min          = 100;
        opt.g_max           = 1e6;                                
  
end

% Derive some options fieldst derived fields
opt = smr_optimize_opt_derive(opt);

%--------------------------------------------------
% PERFORM
%--------------------------------------------------

disp(['Running: ' test_name]);
%
log_fn = fullfile(output_dir, [test_name '_log.txt']);

% Check if performed
if (exist(log_fn, 'file') && ~do_rerun), return; end

% Start log
opt.fid = fopen(log_fn, 'w');

% Perform requested test nr
result = smr_optimize_soma(opt);
save(fullfile(output_dir, [test_name '_result.mat']), 'result');

% Close log
fclose(opt.fid);
end

