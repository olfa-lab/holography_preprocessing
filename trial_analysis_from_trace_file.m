clear all
%% get directory info



trace_dir = 'C:\Users\bnste\Downloads\JG1150\JG1150_trace_files';
trace_file_list = dir(fullfile(trace_dir,'*.mat'));
nfiles = size(trace_file_list,1);


save_file = 'C:\Users\bnste\Downloads\JG1150\JG1150_stat_tables.mat';


% set starting file id in case of appending to previous table
file_id_start = 1;



%% set up file info table

% store file info in this table
fileinfotab = table('Size', [nfiles, 9],'VariableNames',{'recid','mouse','date','filename','onecellstim', 'ncellstim','nframes','ntrials','ngoodtrials'}, 'VariableTypes', {'uint32', 'string','uint32','string', 'uint8','uint32','uint32','uint32','uint32'});
fileinfotab{1:end,'recid'} = uint32(0:nfiles-1)' + file_id_start ;

%% set stat parameters
preCalcPeriod = 25; % must be smaller that omit + fullWindow
postCalcPeriod = 15;

% stim frame is always omitted; these are on either side of that frame
Omitpost = 2;
Omitpre = 2;

%"fullwindow" terms just used for extracting trial windows during reshape,
% not for mean and std calculations
fullWindowPreSize=30;
fullWindowPostSize=30; % 0 will be the first frame post stim


 % set startTrialIdx to control where the trial id numbers start for a
 % given recording
 
 startTrialIdx = 1;

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

save(save_file, 'Fcent_list','prestiminds', 'poststiminds','pattab', 'stattab', 'trialtab', 'filetab', 'preCalcPeriod', 'postCalcPeriod', 'Omitpost', 'Omitpre', 'fullWindowPreSize', 'fullWindowPostSize','maskinds');
