function tif_sequence_pattern = tif_pattern_from_initial(init_name)
% Create a general pattern for Jon's 2p experiment tif files from an
% initial name

[filepath,name,ext] = fileparts(init_name);

splitname = strsplit(name,'_');
num_ind = cellfun(@(x) contains(x, "000"), splitname);
trial_ind = find(num_ind(1:end -1) & num_ind(2:end));
if isempty(trial_ind)
    error(strcat("Can't find pattern for initial video name ", init_name));
end

vid_ind = trial_ind+1;

splitname{vid_ind} = '*';

name = join(splitname, '_');

tif_sequence_pattern = fullfile(filepath, strcat(name,ext));
tif_sequence_pattern = tif_sequence_pattern{:};

end

