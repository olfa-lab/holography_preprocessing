clear all
%% get directory info



trace_dir = 'C:\Users\bnste\Downloads\new_analysis\behavior_end_extracts';
%trace_dir = 'C:\Users\bnste\Downloads\JG1150\stim_extracts_redone_rois';
trace_file_list = dir(fullfile(trace_dir,'JG1150*.mat'));
nfiles = size(trace_file_list,1);

save_file = 'C:\Users\bnste\Downloads\new_analysis\JG1150_behavior_end_tables.mat';
%save_file = 'C:\Users\bnste\Downloads\JG16053\JG16053_tables_nobg.mat';
%save_file = 'C:\Users\bnste\Downloads\JG1150\JG1150_stims_newROIs.mat';


% set starting file id in case of appending to previous table
file_id_start = 1;

% select whether to detrend or not
detrend_data = false;

% set background mask index, if any, and whether you'd like to
% background-subtract
bg_mask_ind = 116; % need to define this if using bg_subtract
bg_subtract = false;
%% set up file info table

% store file info in this table
fileinfotab = table('Size', [nfiles, 9],'VariableNames',{'recid','mouse','date','filename','onecellstim', 'ncellstim','nframes','ntrials','ngoodtrials'}, 'VariableTypes', {'uint32', 'string','uint32','string', 'uint8','uint32','uint32','uint32','uint32'});
fileinfotab{1:end,'recid'} = uint32(0:nfiles-1)' + file_id_start ;

%% set stat parameters
preCalcPeriod = 15; % must be smaller that omit + fullWindow
postCalcPeriod = 3;

% stim frame is always omitted; these are on either side of that frame
Omitpost = 1; %default 1
Omitpre = 2; % default 2

%"fullwindow" terms just used for extracting trial windows during reshape,
% not for mean and std calculations
fullWindowPreSize=45;
fullWindowPostSize=75; % 0 will be the first frame post stim


% recording info
fps=30; % recording speed
tau_rise = 0.5 ; % spike rise time in seconds for fluorescence indicator (used to design filter);


 % set startTrialIdx to control where the trial id numbers start for a
 % given recording
 
 startTrialIdx = 1;
 
 % if a stim pattern has any stim site with a max cell overlap less than
 % this number of pixels, omit that site from the pattern. default in
 % trial_stat_tables.m is 20
 
 minStimOverlap = 50;

%% loop over files, calculate stats, concatenate stats tables, and fill in file info

% store stim-centered fluorescence traces in a Map keyed on recid
Fcent_list = containers.Map('KeyType',class(fileinfotab.recid(1)), 'ValueType','any');



for loopInd=1:nfiles
    disp(sprintf("Extracting stats tables from trace file %d / %d",loopInd ,nfiles));
    trace_file = trace_file_list(loopInd);
    expfile = load(fullfile(trace_file.folder, trace_file.name));
    fileinfotab.filename(loopInd) = trace_file.name;
    recid = fileinfotab.recid(loopInd);
    fileinfotab.date(loopInd) = expfile.date;
    fileinfotab.mouse(loopInd) = expfile.mouse;
    
   

    % create stats tables trialstattab, trialinfotab, patterninfotab and
    % get oneCellStim parameter
    trial_stat_tables;

    Fcent_list(recid) = Fcentered;
    
    fileinfotab.onecellstim(loopInd) = oneCellStim;
    fileinfotab.ncellstim(loopInd) = size(unique(patterninfotab.cell),1);
    fileinfotab.ntrials(loopInd) = ntrialsorig;
    fileinfotab.ngoodtrials(loopInd) = ntrialsgood;
    fileinfotab.nframes(loopInd) = nframes;
    
    patterninfotab{1:end,'recid'} = recid;
    trialstattab{1:end,'recid'} = recid;
    trialinfotab{1:end,'recid'} = recid;
    
    
    if loopInd==1
        patterninfotabfull = patterninfotab;
        trialstattabfull = trialstattab;
        trialinfotabfull = trialinfotab;
    else
        patterninfotabfull = vertcat(patterninfotabfull,patterninfotab);
        trialstattabfull = vertcat(trialstattabfull, trialstattab);
        trialinfotabfull = vertcat(trialinfotabfull, trialinfotab);
    end
    

    % update startTrialIdx for next recording NOTE: giving up on this since
    % we can just join on trial and recind together
    %startTrialIdx = startTrialIdx + size(trialinfotab,1);
end

pattab = patterninfotabfull;
stattab = trialstattabfull;
trialtab = trialinfotabfull;
filetab = fileinfotab;

% save mask info for spatial reconstruction
maskinds = expfile.maskinds;

% use this to include stim centered traces (large, comparitively)
save(save_file, 'detrend_data','bg_subtract','Fcent_list','prestiminds', 'poststiminds','pattab', 'stattab', 'trialtab', 'filetab', 'preCalcPeriod', 'postCalcPeriod', 'Omitpost', 'Omitpre', 'fullWindowPreSize', 'fullWindowPostSize','maskinds');

% use this to exclude them
%save(save_file,'detrend_data','bg_subtract','prestiminds', 'poststiminds','pattab', 'stattab', 'trialtab', 'filetab', 'preCalcPeriod', 'postCalcPeriod', 'Omitpost', 'Omitpre', 'fullWindowPreSize', 'fullWindowPostSize','maskinds');
