function set_me_up(gitDir,fsDir)

setenv('PATH', [fsDir '/bin' getenv('PATH')]);
setenv('FREESURFER_HOME',fsDir)
addpath(genpath(fullfile(fsDir, 'matlab')));

addpath(genpath(fullfile(gitDir, 'cvncode'))); % https://github.com/cvnlab/cvncode
addpath(genpath(fullfile(gitDir, 'knkutils'))); % https://github.com/cvnlab/knkutils
addpath(genpath(fullfile(gitDir, 'GLMdenoise'))); % https://github.com/cvnlab/GLMdenoise
addpath(genpath(fullfile(gitDir, 'nsdcode'))); % https://github.com/cvnlab/nsdcode

end


https://github.com/putiw/glm-Localizer-Tutorial-Code
https://github.com/cvnlab/cvncode
https://github.com/cvnlab/knkutils
https://github.com/cvnlab/GLMdenoise
https://github.com/cvnlab/nsdcode

https://drive.google.com/file/d/1Ot5_QWl6whpB5qEXUEwQuzUJyis-oKJp/view?usp=sharing