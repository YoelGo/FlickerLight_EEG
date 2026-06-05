function [EEGchFixed] = fixChannels(cfg, EEGdat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix (interpolate) bad EEG channels over specified time ranges.
%
% CFG FIELDS:
%   cfg.badch    = cell array of channel names  {'F6', 'Cz', 'F8'}
%   cfg.section  = cell array of time ranges, one entry per channel.
%                  Each entry is either:
%                    'all'  – fix the entire recording
%                    a cell array of [start end] pairs in SECONDS,
%                    where start/end can be a number, 'start', or 'last'
%   cfg.elec     = electrode structure (e.g. standard_1020)
%   cfg.layout   = EEG cap layout
%
% EXAMPLE:
%   cfg.badch   = {'F6',  'Cz',                        'F8'      }
%   cfg.section = {'all', {['start' 120], [350 'last']}, {[120 150]}}
%
%   → F6  : interpolated across the whole recording
%   → Cz  : interpolated from t=0 to t=120 s, and from t=350 s to end
%   → F8  : interpolated from t=120 s to t=150 s
%
% NOTES:
%   • 'start' resolves to the first time-point of the data.
%   • 'last'  resolves to the last  time-point of the data.
%   • For multi-trial (epoched) data every matching trial is repaired;
%     time is matched against each trial's own time axis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(cfg, 'badch') || isempty(cfg.badch)
    warning('fixChannels: no channels specified – returning data unchanged.');
    EEGchFixed = EEGdat;
    return
end
badchannel = cfg.badch;

if any(~ismember(badchannel, EEGdat.label))
    error('fixChannels: one or more channel names not found in data.');
end

if ~isfield(cfg, 'section')
    error('fixChannels: cfg.section is required.');
end
section = cfg.section;

% Allow a single 'all' string to apply to every channel
if ischar(section) && strcmp(section, 'all')
    section = repmat({'all'}, 1, length(badchannel));
end

if length(section) ~= length(badchannel)
    error('fixChannels: cfg.section must have one entry per channel in cfg.badch.');
end

if ~isfield(cfg, 'layout'), error('fixChannels: cfg.layout is required.'); end
if ~isfield(cfg, 'elec'),   error('fixChannels: cfg.elec is required.');   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elec   = cfg.elec;
layout = cfg.layout;

nTrials      = length(EEGdat.trial);
isMultiTrial = nTrials > 1;

% Prepare neighbours once
cfgN         = [];
cfgN.layout  = layout;
cfgN.method  = 'triangulation';
neighbours   = ft_prepare_neighbours(cfgN, EEGdat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper: resolve 'start' / 'last' to numeric seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function t = resolveTime(val, timeAxis)
        if ischar(val)
            switch lower(val)
                case 'start', t = timeAxis(1);
                case 'last',  t = timeAxis(end);
                otherwise
                    error('fixChannels: unknown time keyword "%s". Use "start" or "last".', val);
            end
        else
            t = val;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper: parse one channel's section entry into Nx2 numeric matrix
%%  Returns rows of [tStart tEnd] in seconds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function ranges = parseSection(seg, timeAxis)
        if ischar(seg) && strcmp(seg, 'all')
            ranges = [timeAxis(1), timeAxis(end)];
            return
        end

        % Wrap a bare [start end] pair in a cell so the loop below is uniform
        if isnumeric(seg)
            seg = {seg};
        end

        ranges = zeros(length(seg), 2);
        for k = 1:length(seg)
            pair = seg{k};
            if length(pair) ~= 2
                error('fixChannels: each time range must have exactly 2 elements.');
            end
            ranges(k,1) = resolveTime(pair(1), timeAxis);
            ranges(k,2) = resolveTime(pair(2), timeAxis);
            if isnumeric(ranges(k,1)) && isnumeric(ranges(k,2)) && ranges(k,1) > ranges(k,2)
                error('fixChannels: time range start (%.2f) is after end (%.2f).', ...
                      ranges(k,1), ranges(k,2));
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main loop – one iteration per bad channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(badchannel)

    chan_idx = find(strcmp(badchannel{i}, EEGdat.label));

    % Channel selection: keep all channels except the *other* bad ones
    % (so neighbours of the current channel are still available)
    otherBad = badchannel;
    otherBad(i) = [];

    cfgS         = [];
    cfgS.channel = ft_channelselection([{'all'}, strcat('-', otherBad)], EEGdat.label);

    cfgR             = [];
    cfgR.badchannel  = {badchannel{i}};
    cfgR.neighbours  = neighbours;
    cfgR.method      = 'spline';
    cfgR.elec        = elec;

    % ------------------------------------------------------------------ %
    if isMultiTrial
    % ================================================================== %
    %  MULTI-TRIAL (epoched) data
    %  Time axis is per-trial; ranges are matched against each trial.
    % ================================================================== %
        for tr = 1:nTrials
            trTimeAxis = EEGdat.time{tr};
            ranges     = parseSection(section{i}, trTimeAxis);

            for r = 1:size(ranges, 1)
                tStart = ranges(r,1);
                tEnd   = ranges(r,2);

                % Skip trial if it doesn't overlap with the requested range
                if tEnd < trTimeAxis(1) || tStart > trTimeAxis(end)
                    continue
                end

                cfgS.trials  = tr;
                cfgS.latency = [max(tStart, trTimeAxis(1)), ...
                                min(tEnd,   trTimeAxis(end))];
                tmp = ft_selectdata(cfgS, EEGdat);
                tmp = ft_channelrepair(cfgR, tmp);

                [~, idx1] = min(abs(trTimeAxis - tStart));
                [~, idx2] = min(abs(trTimeAxis - tEnd));
                EEGdat.trial{tr}(chan_idx, idx1:idx2) = tmp.trial{1}(chan_idx, :);
            end
        end

    else
    % ================================================================== %
    %  SINGLE continuous trial
    % ================================================================== %
        timeAxis = EEGdat.time{1};
        ranges   = parseSection(section{i}, timeAxis);

        for r = 1:size(ranges, 1)
            tStart = ranges(r,1);
            tEnd   = ranges(r,2);

            cfgS.latency = [tStart, tEnd];
            tmp = ft_selectdata(cfgS, EEGdat);
            tmp = ft_channelrepair(cfgR, tmp);

            [~, idx1] = min(abs(timeAxis - tStart));
            [~, idx2] = min(abs(timeAxis - tEnd));
            EEGdat.trial{1}(chan_idx, idx1:idx2) = tmp.trial{1}(chan_idx, :);
        end
    end
end

EEGchFixed = EEGdat;
end