function bgmask = addBackgroundMask(maskdir, saveMask)

% use a directory of cell mask image files to generate a random background
% mask that avoids a dilated region of the original masks

% set saveMask to 1 to save the mask to the mask directory with a _bg tag
% saveMask defaults to 0


if nargin < 2
    saveMask = 0;
end

% fraction of total image pixels to use for random selection
fracBg = 1e-2;

% number of pixels to dilate masks for avoidance in bg mask
dilatePix = 10;

maskArray = loadMasks(maskdir);
allMasks = max(maskArray,[],3);
maskArraySize = size(maskArray);
frameSize = maskArraySize(1:2);
nMasks = maskArraySize(3);
framePix = prod(frameSize);
numBg = floor(framePix * fracBg);

se = strel('disk',dilatePix,8);
allMasksDilated = imdilate(allMasks,se);


bgInds = randperm(framePix, numBg);
% remove the pixels that intersect a mask
bgInds = bgInds( allMasksDilated(bgInds) < 1);

bgmask = zeros(frameSize);
bgmask(bgInds) = 1;
bgmask = uint8(255*bgmask);



if saveMask==1
    fname = fullfile(maskdir,sprintf("Mask_1%03d_bg.bmp", nMasks+1));
    imwrite(bgmask,char(fname));
end
end