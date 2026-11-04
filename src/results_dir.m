function folder = results_dir()
%RESULTS_DIR  Folder for generated figures (<project>/results), created if missing.

folder = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results');
if ~exist(folder, 'dir')
    mkdir(folder);
end
end
