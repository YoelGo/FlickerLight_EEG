function ftdat = LSL2ft(EEGdat)

ftdat.trial{1}   = EEGdat.time_series;
ftdat.time{1}    = EEGdat.time_stamps(:)';
ftdat.label      = cellfun(@(c) c.label, ...
                    EEGdat.info.desc.channels.channel, 'UniformOutput', false);
ftdat.fsample    = EEGdat.info.effective_srate;
ftdat.sampleinfo = [1 size(ftdat.trial{1},2)];

ftdat = ft_datatype_raw(ftdat);


end
