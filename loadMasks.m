function [mask_array,spotidx] = loadMasks(path)
%[mask_array,Spotidx] = loadMasks(path)
%   loadMasks() imports .bmp binary masks by file name and outputs a 
%   logical tensor and linear index cell array. 
%
%   loadMasks(path) specifies optional path to folder containing the masks
%   
%   JG 2018

fpathPat = char(path);


mask_list = dir(fullfile(fpathPat,'*.bmp'));
spotidx = cell(length(mask_list),1);
width = 0; height=0;
for idx = 1:length(mask_list)
    temp = imread(fullfile(fpathPat,mask_list(idx).name));
    temp = logical(temp./max(max(temp)));
    if width ==0
        [height, width] = size(temp);
        mask_array = NaN(height,width,length(mask_list));
    end
    mask_array(:,:,idx) = temp;
    spotidx{idx} = find(temp); % *Note* retained Spotidx for legacy use
end
mask_array = logical(mask_array);

end