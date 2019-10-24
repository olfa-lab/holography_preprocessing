function frameidx = get_stim_frames(fpathH5, fnameH5)

%% expt file info


% Getting the signals from the HDF5 file, "M" keeps the data in
% Trials structure and "m" concatenates all trials to one colom.
[ M, m ] = HDF5reader( fpathH5,fnameH5, 0,'frame_triggers','sniff','lick1','lick2');
[ M.packet_sent_time, M.sniff_samples, m.packet_sent_time,...
    m.sniff_samples] = HDF5Eventsreader( fpathH5,fnameH5);

[ m.sniff ] = Check_sniff_data( m );



% This is the local time of the experiment according to sniff recordings
m.time=1:size(m.sniff,1);
m.ftrig_shift=m.frame_triggers(1:end)-single(m.packet_sent_time(1))...
    -single(m.sniff_samples(1));

% Number of trials
m.numTrials = size(m.shutter_onset,1);
m.trial = 1:m.numTrials;

% Timing the shutter onset (shifting the shutter onset time to sniff time)
m.shonset_shift=m.shutter_onset-(m.packet_sent_time(1))...
    -(m.sniff_samples(1));

% Remove stims and trials before the imaging started
m.shutter_timing=zeros(1,size(m.sniff,1));
m.includedTrials = nonzeros(double(m.trial') .* double(m.shonset_shift>0));
m.excludedTrials = setdiff(m.trial,m.includedTrials);
m.shutter_timing(nonzeros(double(m.shonset_shift).*double(m.shonset_shift>0)))=1;

% Find closest frame to trigger onset
onsets=find(m.shutter_timing==1);
diffs=min(abs(m.ftrig_shift-find(m.shutter_timing==1)));
for idx1 = 1:size(nonzeros(m.shonset_shift>0),1)
    frameidx(idx1)=find(abs(m.ftrig_shift-onsets(idx1))==diffs(idx1),1); % This is the closest frame to the shutter onset
end

end
