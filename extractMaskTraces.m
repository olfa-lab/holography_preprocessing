function maskTraces = extractMaskTraces(mov, masks, spotidx)
% extract projected intensity traces from a movie using a mask array
%   mov:  height x width x T
%   masks: height x width x nMasks. expected to be sparse.
%   maskTraces: nMasks x T


if nargin < 3
    findspot = true;
else
    findspot = false;
end


[height,width, T] = size(mov);
nMasks = size(masks,3);

maskTraces = zeros(nMasks,T);
if findspot
    for mInd=1:nMasks
        [i,j] = find( masks(:,:, mInd));
        %maskTraces(mInd,:) = sum(sum(mov(i,j,:),2),1);
        tmp = reshape(mov(i,j,:), numel(i)*numel(j),[]); % pixels x time
        maskTraces(mInd,:) = mean(tmp,1);
        
    end
else
    for mInd=1:nMasks
        [i,j] = ind2sub([height width],spotidx{mInd});
        %maskTraces(mInd,:) = sum(sum(mov(i,j,:),2),1);
        tmp = reshape(mov(i,j,:), numel(i)*numel(j),[]); % pixels x time
        maskTraces(mInd,:) = mean(tmp,1);
    end
end