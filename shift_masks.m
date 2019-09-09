function [] = shift_masks(sourcef, targetf, sourceMaskDir, outMaskDir)



%% Load source template, target templates, and source masks
disp(fprintf("Registering mask source reference %s to mask target reference %s", sourcef, targetf));
disp(fprintf("Source mask dir %s ; output mask dir %s",sourceMaskDir, outMaskDir));

mfiles = dir(fullfile(sourceMaskDir, '*.bmp'));
nMasks = length(mfiles);

source = read_file(sourcef);
source = double(source(:,:,1));
source = log(source+ 1 - min(source(:)));  

target = read_file(targetf);
target = double(target(:,:,1));
target = log(target+1 - min(target(:)));

srcMasks = zeros( [size(source) nMasks]);
for i=1:nMasks
   mf = mfiles(i);
   srcMasks(:,:,i) = imread(fullfile(mf.folder,mf.name));
end

mmax = max(max(srcMasks)); % src mask maxima

%% set up parpool

% create a local cluster object
pc = parcluster('local')

% cluster job handling
if numel(getenv('SLURM_DIR')) > 0
    % explicitly set the JobStorageLocation to the temp directory that was created in your sbatch script
    pc.JobStorageLocation = getenv('SLURM_DIR')

    % start the matlabpool with maximum available workers
    % control how many workers by setting ntasks in your sbatch script
    parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))


%% get shifts


options_nonrigid = NoRMCorreSetParms('d1',size(target,1),'d2',size(target,2),...
    'grid_size',[64,64],'mot_uf',4,'bin_width',50,'max_shift',30,...
    'max_dev',3,'us_fac',50, 'shifts_method','cubic');
tic;

%[M2,shifts2,template2] = normcorre_batch(targetStack,options_nonrigid,source);
%tmpTarget = repmat(targetStack(:,:,i),1,1,nMasks);
[M_final,shifts,template,options,col_shift] = normcorre(source,options_nonrigid,target);

toc


%% Apply shifts to masks

outMasks = zeros(size(srcMasks));
tic;

for mInd=1:nMasks
    outMasks(:,:,mInd) = apply_shifts(srcMasks(:,:,mInd), shifts, options_nonrigid); %,'col_shift',col_shift);
end

outMasks = (outMasks > (.5 * mmax )) .* mmax; % threshold to half mask value to handle interpolation
toc
%% Save output masks
mkdir(outMaskDir);
for mInd=1:nMasks
    mf = mfiles(mInd);
    imwrite(outMasks(:,:,mInd), fullfile(outMaskDir,mf.name))
end

disp(['Registered masks written to output directory ' outMaskDir])

end
