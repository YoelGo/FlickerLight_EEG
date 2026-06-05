%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Flicker Light %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [env] = setupEnviroment11(cfg)

projPath = cfg.project_folder;
dataPath = cfg.data_folder;
ftPath   = cfg.fieldtrip_path;

%%%%% these paths need to be adjusted %%%%%
fieldtrip_path      = '';
layout_path         = [projPath '\extras\resources\LandauLab_BP64ch_Layout.mat'];
analysisdir        = [projPath 'analysis\'];
datadir             = dataPath;

env = [];
%%
%%%%% initiate fieldtrip %%%%%
try
    env.paths.fieldtrip_path = ftPath;
    addpath(ftPath);
    addpath([fieldtrip_path 'external\xdf\']); % add specific LSL (.xdf) file functions
    ft_defaults;
catch
    error('Something is wrong with the fieldtirp path you entered.')
end
%%
%%%%% innitiate env (envelope) variable %%%%%
%%% set paths
env.paths.rawData   = datadir;
env.paths.cleanData = [analysisdir 'preproc\clean\'];
env.paths.artifacts = [analysisdir 'preproc\artifacts\'];

%%% set data variables
%env.data.cleanID         = cellfun(@(x) regexprep(x.ID, '_.*', ''), env.data.dfLAVI, 'UniformOutput', false);
env.data.rawFiles        = dir(fullfile(env.paths.rawData, '**', '*.xdf'));
env.data.rawFileNames    = {env.data.rawFiles.name};    
env.data.cleanFiles      = dir(fullfile(env.paths.cleanData, '**', '*.mat'));
env.data.cleanFileNames  = {env.data.cleanFiles.name};    

%%% set EEG variables
env.EEG.elec        = ft_read_sens([fieldtrip_path 'template\electrode\standard_1020.elc']);
env.EEG.fsample     = 500;

% load EEG layout
cfg          = [];
cfg.layout   = fullfile(layout_path);
env.EEG.lay  = ft_prepare_layout(cfg);


end
