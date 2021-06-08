%
% load_mt.m
%
% Set paths to additional utilities (such as rotation scripts) found in the
% repository surfacevel2strain.
% Get the repo by cloning it at the same location as mtbeach:
% git clone --depth=1 https://github.com/carltape/surfacevel2strain
%
% load_mt.m can be called from a matlab startup.m script,
% so that you never need to execute it.

% path to the directory containing mtbeach and surfacevel2strain
% it is safest to specify an absolute path such as
%repopath = strcat(getenv('HOME'),'/REPOSITORIES/');
% this relative path will only work if you run Matlab from mtbeach/matlab/
repopath = '../../';

addpath(strcat(repopath,'surfacevel2strain/matlab/util_basic'));
addpath(strcat(repopath,'surfacevel2strain/matlab/util_euler'));
addpath(strcat(repopath,'mtbeach/matlab'));
addpath(strcat(repopath,'mtbeach/matlab/Vhat'));
