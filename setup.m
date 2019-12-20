function setup(do_restoredefaultpath)
%
if (nargin < 1), do_restoredefaultpath = 0; end
%
% Restores the default paths and adds all relevant subdirs to the
% path. Run this when you start MATLAB to use the code.
%

subfolders_list = {...
    'estimates', ...
    'figures', ...
    'model', ...
    'optimization', ...
    'optimization/soma', ...
    'tools', ...
    'tools/uvec', ...
   };

t = fileparts(mfilename('fullpath'));

if (do_restoredefaultpath)
    disp('Restoring default path');
    restoredefaultpath;
end
%
for c_package = 1:numel(subfolders_list)
    addpath(fullfile(t, subfolders_list{c_package}), '-end');
end
%
disp (' ')
disp 'Paths are setup.';
