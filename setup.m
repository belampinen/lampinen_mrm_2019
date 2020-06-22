function setup
%
% Restores the default path and adds all relevant subdirs to the
% path. Run this when you start MATLAB to use the code.
%

%%% Setup MDM framework
%
run([pwd '/../md-dmri/setup_paths.m'])

%%% Setup SMR-specific code
%
subfolders_list = {...
    'dependencies', ...
    'estimates', ...
    'figures', ...
    'model', ...
    'optimization', ...
    'optimization/soma', ...
    'tools', ...
    'tools/uvec', ...
   };

t = fileparts(mfilename('fullpath'));

%
for c_package = 1:numel(subfolders_list)
    addpath(fullfile(t, subfolders_list{c_package}), '-end');
end
%
disp (' ')
disp 'Paths are setup.';
