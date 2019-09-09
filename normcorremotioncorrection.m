function normcorremotioncorrection(movname,tempname, redchannel, replace)

% redchannel=0 if there isn't a structural channel to remove, 1 if there is
% replace=1 replace existing output file if =1

if nargin < 3
   redchannel = 0;
   replace = 1;
end
if nargin < 4
   replace = 1;
end


movname
tempname
redchannel

savepart = char(movname);
savepart = string(savepart(1:end-4));
savefile= char(strcat(pwd,'/aligned/',savepart,'_aligned.tif'))

if isfile(savefile) & ~replace
   disp('Output file already exists. Exiting due to "replace" flag value of 0.')
   exit
end

addpath('/gpfs/home/stetlb01/normcorre-matlab/')
% name = ['/experiment/TwoPhoton/2P_Detection/JG8432/170120/Template/' name];
% tempname = '/experiment/TwoPhoton/2P_Detection/JG8432/170120/Template/green/Template_green.tif';

tic; Y = read_file(movname); toc; % read the file (optional, you can also pass the path in the function instead of Y)
Y=Y(:,:,1:(redchannel+1):end);
Y = double(Y);      % convert to double precision 
T = size(Y,ndims(Y));
template = read_file(tempname);

if ndims(template)==3
   template = template(:,:,1);
end



options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),...
    'grid_size',[32,32],'mot_uf',4,'bin_width',50,'max_shift',15,...
    'max_dev',3,'us_fac',50);
tic; 
[M2,shifts2,template2] = normcorre_batch(Y,options_nonrigid,template); 
toc


M2 = cast(M2, 'single');

saveastiff(M2,savefile);


