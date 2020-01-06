%% experiment trace file should be loaded as 'expfile'
% starting trial index recorded as 'startTrialIdx'. this will  be the
% offset for the trial indices


% if ~exist('plotOn')
%    plotOn = false;
% end

frameidx = expfile.frameidx;
F = expfile.F;


if ~exist('startTrialIdx', 'var')
    startTrialIdx = 1;
end
if ~exist('minStimOverlap', 'var')
    minStimOverlap = 20; %minimum required pixel overlap
end
if ~exist('detrend_data','var')
    detrend_data = false;
end
if ~exist('bg_subtract','var')
    bg_subtract = false;
end

% match stim sites in SpotMat, Spotidx (and used in patternTrials) to
% masks
maskIdx=zeros(size(expfile.Spotidx));
maskMatchMeas = zeros(size(expfile.Spotidx));
nOverlaps = {};
for i=1:numel(maskIdx)
    thisSpot = expfile.SpotMat{i}; % this is the stim spot
    overlaps = cellfun(@(x) sum(thisSpot(x)), expfile.maskinds); %expfile.maskinds is a cell array of the  linear index lists for the cell masks
    overlaps(bg_mask_ind) = -1; % prevent matching with bg mask
    nOverlaps{i} = overlaps;
    [maskMatchMeas(i) , maskIdx(i)] = max(overlaps);
end
maskIdx = uint16(maskIdx);



% % need to handle non-unique matches - make sure we get unique maskIdx. (is this always true?
% % it turns out no, there are repeats sometimes)
% [uniqueM i j] = unique(maskIdx,'first');
% indexToDupes = find(not(ismember(1:numel(maskIdx),i)));
% tryCount=0;
% if numel(indexToDupes) > 0
%     dupe_vals = unique(maskIdx(indexToDupes));
%     disp(sprintf("Warning: overlapping stim sites / masks in mouse %s, date %s", expfile.mouse, expfile.date));
%     for z=1:numel(dupe_vals)
%        disp(sprintf("Multiple mask assignments to mask %d",dupe_vals(z)));
%     end
% end
% while numel(indexToDupes) > 0 && tryCount < 1000
%     dupe_vals = unique(maskIdx(indexToDupes));
%     for dval_ind=1:numel(dupe_vals)
%         dval = dupe_vals(dval_ind);
%         val_inds = find(maskIdx == dval);
%         [~, bestIdx] = max(maskMatchMeas(val_inds));
%         val_inds(bestIdx) = [];
%         avail_inds = find(~ismember(1:numel(nOverlaps{1}),maskIdx));
%         for kind=1:numel(val_inds)
%             k = val_inds(kind);
%             laps = nOverlaps{k};
%             [v,vinds] = sort(laps, 'descend');
%             vinds = vinds(ismember(vinds,avail_inds)); v = v(ismember(vinds,avail_inds));
%             maskIdx(k) = vinds(1);
%             maskMatchMeas(k) = v(1);
%         end
%     end
%     [uniqueM i j] = unique(maskIdx,'first');
%     indexToDupes = find(not(ismember(1:numel(maskIdx),i)));
%     tryCount = tryCount + 1;
% end
% if numel(indexToDupes) > 0
%     error("Cannot create 1-1 matching of stim sites and masks!");
% end


% create vector recording the stim mask number for each
% trial

if isfield(expfile,'patternTrials')
    oneCellStim = true;
else
    oneCellStim = false;
end

if oneCellStim
    [patternTrialsUnroll,~] = find(cell2mat(cellfun(@(x) ismember(1:size(frameidx,2),x), expfile.patternTrials, 'UniformOutput',  false)'));
%     maskTrials = maskIdx(patternTrialsUnroll);
else
    patternTrialsUnroll = ones(1, numel(frameidx));
%     maskTrials = repmat(maskIdx,[1 numel(frameidx)]);
end
    patternTrialsUnroll = uint16(patternTrialsUnroll);
%     maskTrials = uint16(maskTrials);


%% eliminate trials overlapping start/end, and interpolate through stim periods

ntrialsorig = size(frameidx,2);
nframes = size(F,1);

badframes = reshape( (-Omitpre:Omitpost)' + frameidx,1,[] ); % frames within the bad period surrounding a stim
badframes = max(badframes,1); badframes = min(badframes,nframes);
fullwindow = (-fullWindowPreSize:fullWindowPostSize)' + frameidx;

% eliminate bad trials (window overlaps start/end of experiment)
movlen = size(F,1);
badtrials = any( (fullwindow <= 0) | (fullwindow > movlen), 1);
frameidx(badtrials) = [];
patternTrialsUnroll(badtrials) = [];

ntrialsgood = size(frameidx,2);

% linearly interpolate through stim frames
F(badframes,:) = nan;
F = fillmissing(F,'linear',1);

%% detrend data, if desired

if detrend_data
    for icell=1:size(F,2)
        F(:,icell) = detrend(F(:,icell),2); % quadratic detrend
    end
end




%% reshape F to trial-centered form

fullwindow = [(-fullWindowPreSize:-1) (0:fullWindowPostSize)]' + frameidx; % redo window creation
fullwindow = reshape(fullwindow,1,[]);
fullwindow = min(fullwindow,size(F,1));

Fcentered = reshape(F(fullwindow,:),fullWindowPreSize + fullWindowPostSize + 1,[],size(F,2));

% get indices of pre-post periods in the new stim-centered array
prestiminds = fullWindowPreSize + (-preCalcPeriod+1:0) - Omitpre;
poststiminds = fullWindowPreSize + 1 + Omitpost + (1:postCalcPeriod);


%% background subtract, if desired
% NOTE: we fit on the pre-trial periods to avoid correlations due to
% network activity during spiking

if bg_subtract
    
    if ~exist('bg_mask_ind','var')
        error("Variable 'bg_mask_ind' must be defined if bg_subtract is used")
    end
    non_bg_inds = setdiff(1:size(Fcentered,3),bg_mask_ind);
    
    trialFitInds = 1:size(Fcentered,1); % full trial fit
    %trialFitInds = 1:fullWindowPreSize-Omitpre % fit solely to pretrial
    fsize = size(Fcentered);
    Fbg = reshape(Fcentered(trialFitInds,:,bg_mask_ind),[],1);  % want to fit BEFORE stim, when traces aren't as correlated
    
    % lowpass filter on background signal
    bgfilt = designfilt('lowpassfir','FilterOrder',8,'PassbandFrequency',0.1/tau_rise,'StopbandFrequency',0.2/tau_rise,'SampleRate',fps);
    Fbg = filtfilt(bgfilt,Fbg);
    
    Fnobg = reshape(Fcentered(trialFitInds,:,non_bg_inds),[],numel(non_bg_inds));
    Fbgstd = std(Fbg); Fnobgstd = std(Fnobg,0,1);
    Fbgmean = mean(Fbg); Fnobgmean = mean(Fnobg,1);
    Fbg = (Fbg - Fbgmean) ./ Fbgstd; Fnobg = (Fnobg - Fnobgmean) ./ Fnobgstd;
    
    % least squares
    beta = zeros(1,numel(non_bg_inds));
    for j=1:size(beta,2)
        beta(:,j) = Fbg \ Fnobg(:,j);
    end

%     % nonneg least squares
%     betann = zeros(1,numel(non_bg_inds));
%     for j=1:size(betann,2)
%         betann(:,j) = lsqnonneg(Fbg,Fnobg(:,j ));
%     end
    
    
    
    % convex optimization
    %betaopt = fmincon(@(beta) overshootPenaltyFunc(beta,Fbg,Fnobg),betann,[],[],[],[],0,Inf);
    
    % apply linear model to whole traces
    Fbg = reshape(Fcentered(:,:,bg_mask_ind),[],1);
    Fbg = filtfilt(bgfilt,Fbg);
    %Fbg_list(recid) = Fbg;
    Fbg = (Fbg - Fbgmean) ./ Fbgstd;
    Fnobg = reshape(Fcentered(:,:,non_bg_inds),[],numel(non_bg_inds));
    Fnobg = (Fnobg - Fnobgmean) ./ Fnobgstd;
    Fpred = Fbg*beta;

    Fresid = ones(size(Fpred,1),size(Fcentered,3));
    Fresid(:,non_bg_inds) = Fnobg - Fpred;
    Fresid(:,non_bg_inds) = (Fresid(:,non_bg_inds) .* Fnobgstd) + Fnobgmean;
    Fresid(:,bg_mask_ind) = Fbgmean;
    
    Fcentered = reshape(Fresid, fsize);
end



%% get stim-centered stats (dff, response snr, response type,..)
prestimmean = single(mean(reshape(Fcentered(prestiminds,:,:),numel(prestiminds),[],size(F,2)),1));
poststimmean = single(mean(reshape(Fcentered(poststiminds,:,:),numel(poststiminds),[],size(F,2)),1));
prestimstd = single(std(reshape(Fcentered(prestiminds,:,:),numel(prestiminds),[],size(F,2)),1));
poststimstd = single(std(reshape(Fcentered(poststiminds,:,:),numel(poststiminds),[],size(F,2)),1));
tscore = single(  (poststimmean - prestimmean) ./ sqrt( (prestimstd.^2 ./ size(prestiminds,2)) + (poststimstd.^2 ./ size(poststiminds,2)) ));


dff = single(squeeze((poststimmean - prestimmean) ./ prestimmean));
respsnr = single(squeeze( (poststimmean - prestimmean) ./ prestimstd)); % dF / stddev of prestim measurement period
resptype = int8(reshape( discretize(respsnr(:), [-inf -1 1 inf]) - 2, size(respsnr))); % -1 for inhibition, 1 for response, 0 otherwise
dffmean = squeeze(mean(dff,1));

% get background dff response
spontpreinds = 1:floor(fullWindowPreSize/2);
spontpostinds = (floor(fullWindowPreSize/2)+1) : (fullWindowPreSize - 1);
spontpremean = single(mean(reshape(Fcentered(spontpreinds,:,:),numel(spontpreinds),[],size(F,2)),1));
spontpostmean = single(mean(reshape(Fcentered(spontpostinds,:,:),numel(spontpostinds),[],size(F,2)),1));
spontprestd = single(std(reshape(Fcentered(spontpreinds,:,:),numel(spontpreinds),[],size(F,2)),1));
spontpoststd = single(std(reshape(Fcentered(spontpostinds,:,:),numel(spontpostinds),[],size(F,2)),1));
spontdff = single(squeeze((spontpostmean - spontpremean) ./ spontpremean));
spontrespsnr = single(squeeze( (spontpostmean - spontpremean) ./ spontprestd)); % dF / stddev of prestim measurement period
spontresptype = int8(reshape( discretize(spontrespsnr(:), [-inf -1 1 inf]) - 2, size(spontrespsnr))); % -1 for inhibition, 1 for response, 0 otherwise

%% get random dff response
randframes = sort(randsample(fullWindowPreSize+1:size(F,1)-fullWindowPostSize-1, numel(frameidx), true));

fullwindow = [(-fullWindowPreSize:-1) (0:fullWindowPostSize)]' + randframes; % redo window creation
fullwindow = reshape(fullwindow,1,[]);
fullwindow = min(fullwindow,size(F,1));

Fcentrand = reshape(F(fullwindow,:),fullWindowPreSize + fullWindowPostSize + 1,[],size(F,2));

% get indices of pre-post periods in the new stim-centered array
randpreinds = fullWindowPreSize + (-preCalcPeriod+1:0) - Omitpre;
randpostinds = fullWindowPreSize + 1 + Omitpost + (1:postCalcPeriod);

randpreinds = 1:floor(fullWindowPreSize/2);
randpostinds = (floor(fullWindowPreSize/2)+1) : (fullWindowPreSize - 1);
randpremean = single(mean(reshape(Fcentrand(randpreinds,:,:),numel(randpreinds),[],size(F,2)),1));
randpostmean = single(mean(reshape(Fcentrand(randpostinds,:,:),numel(randpostinds),[],size(F,2)),1));
randprestd = single(std(reshape(Fcentrand(randpreinds,:,:),numel(randpreinds),[],size(F,2)),1));
randpoststd = single(std(reshape(Fcentrand(randpostinds,:,:),numel(randpostinds),[],size(F,2)),1));
randdff = single(squeeze((randpostmean - randpremean) ./ randpremean));
randrespsnr = single(squeeze( (randpostmean - randpremean) ./ randprestd)); % dF / stddev of prestim measurement period
randresptype = int8(reshape( discretize(randrespsnr(:), [-inf -1 1 inf]) - 2, size(randrespsnr))); % -1 for inhibition, 1 for response, 0 otherwise
%% create stat and trial info tables

[triallabel, celllabel] = ind2sub(size(dff), 1:numel(dff));

trialstattab = table( -1 + startTrialIdx + uint32(triallabel(:)), uint16(celllabel(:)),  dff(:), respsnr(:), resptype(:), prestimmean(:), poststimmean(:), prestimstd(:), poststimstd(:), tscore(:), spontdff(:), spontrespsnr(:), spontresptype(:),randdff(:), randrespsnr(:), randresptype(:),'VariableNames', {'trial','cell','dff','respsnr','resptype', 'premean','postmean','prestd','poststd', 'tscore','spontdff','spontrespsnr','spontresptype','randdff','randrespsnr','randresptype'} ); %, 'VariableTypes',{'uint16','uint16','single','single','int8'});
trialinfotab = table( -1 + startTrialIdx + uint32((1:size(dff,1))'), patternTrialsUnroll(:), uint32(frameidx(:)),'VariableNames', {'trial','pattern', 'frame'});
if oneCellStim
    patterninfotab = table( uint16(1:numel(maskIdx))', uint16(maskIdx)', uint32(maskMatchMeas(:)), 'VariableNames',{'pattern','cell', 'matchpix'});
else
    patterninfotab = table( uint16(ones(numel(maskIdx),1)), uint16(maskIdx(:)) , uint32(maskMatchMeas(:)), 'VariableNames',{'pattern','cell', 'matchpix'});
end

patterninfotab = patterninfotab(patterninfotab.matchpix >= minStimOverlap,:);




