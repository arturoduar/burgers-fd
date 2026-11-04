%BURGERS_SETUP  Put the project's source folders on the MATLAB path.
%   Run once per session from the project root:  burgers_setup
root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, 'src')));
fprintf('burgers_fd: source folders added to the path.\n');
