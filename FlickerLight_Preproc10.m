%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Flicker Light %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% loading in LSL file and preprocessing the data
clc; clear all; close all;
% set relevant paths %
% gitPath = 'C:\Users\yoelgo\Documents\GitHub\FlickerLight_EEG'; % set path of the local GitHub
gitPath = 'C:\Users\ayelet.landau\Dropbox\Analysis\FlickerLight_EEG-main\';
addpath(genpath(gitPath)); % add all subfolders in this project

% a custom function which sets all the paths and prepares the variable 
% enviroment (env) for this experiment
cfg = [];
cfg.data_folder    = 'C:\Users\ayelet.landau\Dropbox\Analysis\FlickerLight_EEG-main\'; % insert the path where your data is
cfg.fieldtrip_path = 'C:\Users\ayelet.landau\Dropbox\Analysis\fieldtrip-20240110\fieldtrip-20240110'; % insert the path of your fieldtrip.

% had to add this one:

cfg.project_folder = gitPath; % keep constant
env = setupEnviroment11(cfg); 

% addpath('C:\Users\ayelet.landau\Dropbox\Analysis\fieldtrip-20240110\fieldtrip-20240110\external\xdf\')

clear gitPath
%% Load single participant
% ***find the data number in the raw files list in env.data.rawFiles***
n = 1;

% find the participant ID
ID = extractBefore(env.data.rawFiles(n).name, '.xdf');
env.ID = ID;
env.paths.curData = [env.data.rawFiles(n).folder '\'];

% load data
filePath = [env.paths.curData  env.data.rawFiles(n).name];
LSLdat = load_xdf(filePath);
% extract the EEG channel and convert to fieldtirp.
ftEEG  = LSL2ft(LSLdat{cellfun(@(x) strcmp(x.info.name, 'actiCHamp-21020490'), LSLdat)})
%% divide into the different LSL streams
% This part depends on the triggers you have in your experiment.
%% find markers with the photodiode
% photodiode signal
cfg = [];
cfg.channel = {'AUX_4'}; % I think it's better to use this channel for segmentation in your case
cfg.demean = 'yes';
photoD = ft_preprocessing(cfg,ftEEG);

fs = photoD.fsample;

% threshold
thr = 0.5 * max(pd.trial{1});

% logical vector: 1 = above threshold
above = pd.trial{1} > thr;

% crossings
onsets  = find(diff([0 above]) == 1);   % crosses from below to above
offsets = find(diff([above 0]) == -1);  % crosses from above to below

trlTable = table(onsets', offsets', zeros(length(offsets),2), ...
    ((offsets-onsets)/fs)',[0 (onsets(2:end) - offsets(1:end-1))/fs]', ...
    'VariableNames',{'beg_sample', 'end_sample', 'offset', 'len', 'Tdif'});

% find block type order
% Here I load the CSV with the block conditions in my pilot, so this part
% is probably not relevant in your case
% csvFile = env.data.csvFiles(contains(string({env.data.csvFiles.name}),ID));
% T = readtable(fullfile(csvFile.folder,csvFile.name));
% 
% valid = ~ismissing(string(T.condition)) & string(T.condition) ~= "";
% [~,idx] = unique(T.block_number(valid),'stable');
% 
% conditions = string(T.condition(valid));
% conditions = conditions(idx);
% 
% % conditions = vector of block conditions, in order
% breaks = trlTable.Tdif(:) > 3.5;
% blockIdx = cumsum(breaks);
% blockIdx(blockIdx == 0) = 1;
% 
% trlTable.condition = conditions(blockIdx);
% 
clear T offsets onsets b idx thr above

% divide into blocks
cfg = [];
cfg.trl = trlTable{:,1:3};
trlDat = ft_redefinetrial(cfg,ftEEG);

%% basic preproc
% first keep the EOG to add later
cfg = [];
cfg.channel = {'AUX_1', 'AUX_2'};
cfg.demean = 'yes';
EOG = ft_preprocessing(cfg,ftEEG);

clc;
cfg            = [];
cfg.channel    = {'all', '-AUX_1', '-AUX_2', '-AUX_3', '-AUX_4'}; % remove EOG and ECG for refferencing
cfg.detrend    = 'yes';
cfg.demean     = 'yes';
cfg.reref      = 'yes';
cfg.refchannel = {'all'}; % A1 and A2 are the ear-clips
% Notch filters to remove 50Hz and harmonics
%cfg.lpfreq      = 70;
%cfg.lpfilter    = 'yes';
cfg.hpfilter    = 'yes';
cfg.hpfreq      = 0.5;
cfg.bsfilter    = 'yes';
cfg.bsfreq      = [49 51; 99 101];
cfg.bsfilttype  = 'but';
cfg.bsfiltord   = 3; % 3rd order
cfg.bsfiltdir   = 'twopass'; % zero-phase

pEEG = ft_preprocessing(cfg, ftEEG); %p(preprocessed)EEG

pEEG.label(1:64) = env.EEG.lay.label(1:64);
% add EoG back
pEEG = ft_appenddata([],pEEG,EOG);


clear Eog

%% Manual: View Data
cfg = [];
cfg.ylim  = [-50 50];
cfg.continuous = 'yes';
cfg.blocksize = 30;

man_artifact = ft_databrowser(cfg,pEEG)

save([env.paths.artifacts  ID '_man_artifact'], "man_artifact"); %%% this line needs to be fixed
save(['C:\Users\ayelet.landau\Dropbox\Analysis\FlickerLight_EEG-main\preproc\artifacts' ID '_man_artifact'], "man_artifact")

%% Manual: Remove Artifacts
cfg = []; 
cfg.artfctdef.reject           = 'nan';
cfg.artfctdef.visual.artifact = man_artifact.artfctdef.visual.artifact;
mEEG = ft_rejectartifact(cfg,pEEG);

cfg = [];
cfg.channel = {'all', '-AUX_1', '-AUX_2', '-AUX_3', '-AUX_4'}; %note, we might need AUX_3 AUX_4
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
%% verify that this last step works later!!!

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
%% this last step needs to be checked!

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

cfg = [];
dat_after_ICA = ft_appenddata(cfg, dat_after_ICA, photoD)

%% save data after preproc
save([env.paths.cleanData ID 'EEG' '_clean.mat'], "dat_after_ICA");
disp(['Saved All Data!']);