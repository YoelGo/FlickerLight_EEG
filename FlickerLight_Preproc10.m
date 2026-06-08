%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Flicker Light %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% loading in LSL file and preprocessing the data
clc; clear all; close all;
% set relevant paths %
gitPath = 'C:\Users\yoelgo\Documents\GitHub\FlickerLight_EEG'; % set path of the local GitHub
addpath(genpath(gitPath)); % add all subfolders in this project

% a custom function which sets all the paths and prepares the variable 
% enviroment (env) for this experiment
cfg = [];
cfg.data_folder    = ''; % insert the path where your data is
cfg.fieldtrip_path = ''; % insert the path of your fieldtrip.

cfg.project_folder = gitPath; % keep constant
env = setupEnviroment11(cfg); 

clear gitPath
%% Load single participant
% ***find the data number in the raw files list in env.data.rawFiles***
n = 4;

% find the participant ID
ID = extractBefore(env.data.rawFiles(n).name, '.xdf');
env.ID = ID;
env.paths.curData = [env.data.rawFiles(n).folder '\'];

% load data
dat1 = load_xdf([env.paths.curData  env.data.rawFiles(n).name]);

% find the EEG stream and convert to fieldtrip format
EEG    = dat1{cellfun(@(x) strcmp(x.info.name, 'actiCHamp-21020490'), dat1)};
ftEEG  = LSL2ft(EEG);
% or use xdf2fieldtrip('file path');
%% divide into the different LSL streams
% This part depends on the triggers you have in your experiment.

%% basic preproc
% first keep the EOG to add later
cfg = [];
cfg.channel = {'AUX_1', 'AUX_2'};
cfg.demean = 'yes';
EOG = ft_preprocessing(cfg,ftEEG);

clc;
cfg            = [];
cfg.channel    = {'EEG', '-AUX_1', '-AUX_2', '-AUX_3', '-AUX_4'} ; % remove EOG and ECG for refferencing
cfg.detrend    = 'yes';
cfg.demean     = 'yes';
cfg.reref      = 'yes';
cfg.refchannel = {'all'}; % A1 and A2 are the ear-clips
% Notch filters to remove 50Hz and harmonics
cfg.lpfreq      = 70;
cfg.lpfilter    = 'yes';
cfg.hpfilter    = 'yes';
cfg.hpfreq      = 0.5;
cfg.bsfilter    = 'yes';
cfg.bsfreq      = [49 51; 99 101];
cfg.bsfilttype  = 'but';
cfg.bsfiltord   = 3; % 3rd order
cfg.bsfiltdir   = 'twopass'; % zero-phase

pEEG = ft_preprocessing(cfg, ftEEG); %p(preprocessed)EEG


% add EoG back
pEEG = ft_appenddata([],pEEG,EOG);

clear Eog

%% Manual: View Data
cfg = [];
cfg.ylim  = [-50 50];
cfg.continuous = 'yes';
cfg.blocksize = 30;

man_artifact = ft_databrowser(cfg,pEEG)

save([env.paths.artifacts  ID '_man_artifact'], "man_artifact");

%% Manual: Remove Artifacts
cfg = []; 
cfg.artfctdef.reject           = 'nan';
cfg.artfctdef.visual.artifact = man_artifact.artfctdef.visual.artifact;
mEEG = ft_rejectartifact(cfg,pEEG);

cfg = [];
cfg.channel = {'all', '-AUX_1', '-AUX_2', '-AUX_3', '-AUX_4'};
mEEG = ft_preprocessing(cfg,mEEG);



%% check for further noisy channels
cfg        = [];
cfg.metric = 'var';  % use by default zvalue method
cfg.method = 'summary'; % use by default summary method
mEEG       = ft_rejectvisual(cfg,mEEG);

%% find channel label from the number (since ft_rejectvisual only gives you channel number)
mEEG.label(62)         % dicsover the channel from the number
%% inspct the bad channel
cfg = [];
cfg.ylim  = [-50 50];
cfg.continuous = 'yes';
cfg.blocksize = 30;
ft_databrowser(cfg,mEEG)

%% remove noisy channels and keep them to be interpolated later.
badCh = {'AF4'}; % Insert the noisy channel name here.

if exist('badCh','var') && ~isempty(badCh{1})
    cfg = [];
    cfg.channel = badCh;
    badChDat = ft_selectdata(cfg,mEEG);

    cfg = [];
    cfg.channel = [{'all'}, strcat('-', badCh)];
    mEEG = ft_selectdata(cfg, mEEG);
else
    warning('No bad channels to interpolate?');
end
%% Run ICA
cfg = [];
cfg.method  = 'runica';
cfg.channel = {'all'}; % EEG no EoG, no badch
cfg.numcomponent = 20;
%cfg.runica.maxsteps = 100;     % only uncomment if the ICA is too long
comp = ft_componentanalysis(cfg, mEEG);
%% view ICA components
% view time seriers and topopraphy of ICs
cfg = [];
cfg.viewmode = 'component';
cfg.allowoverlap = 'yes';
cfg.continuous = 'yes';
cfg.blocksize = 30;
cfg.layout = env.EEG.lay;
ft_databrowser(cfg,comp);

%% reject components
cfg = [];
cfg.component = [1 2]; % insert the component numbers to be removed
dat_after_ICA = ft_rejectcomponent(cfg, comp);

%% reintroduce the bad channels and interpolate them
% reintroduce the bad channels
dat_after_ICA = ft_appenddata([], dat_after_ICA, badChDat);

% interpolate the bad channels using a custom function. see notation iside.
cfg = [];
cfg.badch = badCh;
cfg.section = {'all'};
cfg.elec = env.EEG.elec;
cfg.layout = env.EEG.lay;
cfg.blocksize = 30; 
dat_after_ICA = fixChannels12(cfg, dat_after_ICA);
%% check again the channel noise after interpolation
cfg        = [];
cfg.metric = 'zvalue';  % use by default zvalue method
cfg.method = 'summary'; % use by default summary method
dat_after_ICA = ft_rejectvisual(cfg,dat_after_ICA);

%% view the data again - one last look at the data.
cfg = [];
cfg.ylim  = [-30 30];
cfg.continuous = 'yes';
cfg.blocksize = 30;
man_artifact = ft_databrowser(cfg,dat_after_ICA)


%% save data after preproc
save([env.paths.cleanData ID 'EEG' '_clean.mat'], "dat_after_ICA");
disp(['Saved All Data!']);