function set_me_up(gitDir,fsDir,bidsDir)

setenv('PATH', [fsDir '/bin' getenv('PATH')]);
setenv('FREESURFER_HOME',fsDir)
setenv('SUBJECTS_DIR',[bidsDir '/derivatives/freesurfer'])

addpath(genpath(fullfile(fsDir, 'matlab')));

addpath(genpath(fullfile(gitDir, 'cvncode'))); % https://github.com/cvnlab/cvncode
addpath(genpath(fullfile(gitDir, 'knkutils'))); % https://github.com/cvnlab/knkutils
addpath(genpath(fullfile(gitDir, 'GLMdenoise'))); % https://github.com/cvnlab/GLMdenoise
addpath(genpath(fullfile(gitDir, 'nsdcode'))); % https://github.com/cvnlab/nsdcode

end
