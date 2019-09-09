%% test individual cells

mov_dir = '/gpfs/data/shohamlab/shared_data/jon_2p_data/JG1150/190708/aligned';
mov_pattern = 'JG1150_190708_field2_stim_00001_*_aligned.tif';
mask_dir = '/gpfs/data/shohamlab/shared_data/jon_2p_data/JG1150/190708/masks';
expt_info_file = '/gpfs/data/shohamlab/shared_data/jon_2p_data/JG1150/190708/1150_2_01_D2019_7_8T18_38_59.h5';
pattern_path = '/gpfs/data/shohamlab/shared_data/jon_2p_data/JG1150/190708/JG1150_190708_field2_1_individual_cells';
red_channel = 0;

tic
traces_from_files(mov_dir, mov_pattern, pattern_path, mask_dir, expt_info_file,red_channel);
toc

%% test group stim

mov_dir = '/gpfs/data/shohamlab/shared_data/jon_2p_data/JG1150/190712/aligned';
mov_pattern = 'JG1150_190712_field1_beh_00001_*_aligned.tif';
mask_dir = '/gpfs/data/shohamlab/shared_data/jon_2p_data/JG1150/190712/masks';
expt_info_file = '/gpfs/data/shohamlab/shared_data/jon_2p_data/JG1150/190712/1150_5_01_D2019_7_12T17_48_23_beh.h5';
pattern_path = '/gpfs/data/shohamlab/shared_data/jon_2p_data/JG1150/190712/JG1150_190712_field1_revised_1.mat';
red_channel = 0;

tic
traces_from_files(mov_dir, mov_pattern, pattern_path, mask_dir, expt_info_file,red_channel);
toc