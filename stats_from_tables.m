%% load in experiment stats and info tables

clear all;

% load in all derived tables from this script, if desired
%multi_files = {'C:/Users/bnste/Downloads/JG1150/JG1150_specificity_withrand_workspace.mat' 'C:/Users/bnste/Downloads/JG16271/JG16271_specificity_withrand_workspace.mat'};
multi_files = {'C:\Users\bnste\Downloads\new_analysis\JG1150_behavior_end_workspace.mat' 'C:\Users\bnste\Downloads\new_analysis\JG16271_behavior_end_workspace.mat'};
nmultifiles = numel(multi_files);

%tabfile = 'C:\Users\bnste\Downloads\JG16271\JG16271_specificity_stats_withrand_final.mat';
%tabfile = 'C:\Users\bnste\Downloads\JG16271\JG16271_specificity_stats_3framepost.mat';
%tabfile = 'C:\Users\bnste\Downloads\JG1150\JG1150_specificity_stats_withrand_final.mat';
%tabfile = 'C:\Users\bnste\Downloads\JG16053\JG16053_tables_nobg.mat';
tabfile = 'C:\Users\bnste\Downloads\new_analysis\JG1150_behavior_end_tables.mat';
load(tabfile);

% set recording type (0 for both, 1 for one-cell stim, 2 for multi-cell
% stim)
rec_type = 0;



% movie size
height = 512;
width = 512;
pix_um_ratio = 512/365;

% load cell type indicator list here
% -1 is background, 0 is non-red, 1 is red
celltype = zeros(size(maskinds));
celltype(1:52) = 1; celltype(63) = 1;  % JG1150
%celltype(1:53)=1; %JG16503
%celltype(1:87)=1; %JG16271;

% bg mask index
bg_mask_ind = max(stattab.cell); %JG1150
%bg_mask_ind = nan; %JG16053
%bg_mask_ind =  max(stattab.cell); %JG16271

non_bg_inds = setdiff(1:numel(maskinds),bg_mask_ind);
if ~isnan(bg_mask_ind)
    celltype(bg_mask_ind) = -1;
end

% responders
isresponder = (1:numel(celltype))==0;

isresponder(1:19) = true; % JG1150
%isresponder([6 7 14 15 17 18 20 25 32 34 49]) = true; %JG16053
%isresponder(1:27) = true; %JG16271

responders = find(isresponder);



% get one-cell-stim recids
onestim_recids = filetab.recid(filetab.onecellstim==1);
multistim_recids = filetab.recid(filetab.onecellstim==0);
red_cells = find(celltype == 1);



% filter
switch rec_type
    case 1
        valid_recids = onestim_recids;
    case 2
        valid_recids = multistim_recids;
    otherwise
        valid_recids = filetab.recid;
end

filetab = filetab(ismember(filetab.recid, valid_recids),:);
pattab = pattab(ismember(pattab.recid, valid_recids),:);
trialtab = trialtab(ismember(trialtab.recid, valid_recids),:);
stattab = stattab(ismember(stattab.recid, valid_recids),:);

% colormap for plots (blue->white->red)
nmap = 128;
cmin = .2;
mymap = vertcat([linspace(cmin,1,nmap)' linspace(cmin,1,nmap)' ones(nmap,1)], [ones(nmap,1) linspace(1,cmin,nmap)' linspace(1,cmin,nmap)']); % blue->white->red


%% calculate some useful quantities for later

% get stimmed cells
stimcells = table(unique(pattab.cell), 'VariableNames',{'cell'});


% get table of stimmed patterns with the cells in those patterns for each
% recid
patcellmatch =removevars(pattab, {'matchpix'});

% make sure we only have given patterns (some may have been removed due to
% bad matching
trialtab = outerjoin(trialtab,pattab,'RightVariables',{'pattern'});
trialtab = trialtab(trialtab.pattern_pattab > 0,{'pattern_trialtab','recid','trial','frame'});
trialtab.Properties.VariableNames{'pattern_trialtab'} = 'pattern';

stattab = outerjoin(stattab, trialtab,'Type','Left','LeftKeys',{'trial','recid'}, 'RightKeys',{'trial','recid'},'RightVariables',{'pattern'});
stattab = stattab(stattab.pattern > 0,:);


% % add cellpattern column to stats table (pattern corresponding to cell)
% stattab = outerjoin( stattab ,patcellmatch, 'RightVariables',{'pattern'});
% stattab.Properties.VariableNames{'pattern'} = 'cellpattern';


% add cellstim boolean (whether this cell is stimmed)
stattab = outerjoin(stattab, patcellmatch, 'Type','left', 'RightVariables',{'cell','pattern'}); %  = cell if pattern includes cell for that trial, otherwise 0
stattab.cellstim = (stattab.cell_stattab == stattab.cell_patcellmatch);
stattab = removevars(stattab,{'cell_patcellmatch', 'pattern_patcellmatch'});
stattab.Properties.VariableNames{'cell_stattab'} = 'cell';
stattab.Properties.VariableNames{ 'pattern_stattab'} =  'pattern';
% add celltype column
stattab.celltype = celltype(stattab.cell);

stattab.onestim = ismember(stattab.recid, onestim_recids);


nstimcells = numel(stimcells.cell);
ncells = numel(celltype);



% distance mapping
% get mean x,y position for each mask
xavg = zeros(1,numel(maskinds));
yavg = zeros(1,numel(maskinds));
for i=1:numel(maskinds)
    [yinds, xinds] = ind2sub([height width], maskinds{i});
    xavg(i) = mean(xinds(:));
    yavg(i) = mean(yinds(:));
end


xdistmat = xavg - xavg';
ydistmat = yavg -yavg';
distmat = sqrt(xdistmat.^2 + ydistmat.^2);
[icell, istimcell] = ind2sub([numel(maskinds) numel(maskinds)],1:numel(maskinds)^2);
disttab = table(icell', istimcell', distmat(1:numel(maskinds)^2)',xdistmat(1:numel(maskinds)^2)', ydistmat(1:numel(maskinds)^2)' ,'VariableNames',{'viewcell','stimcell','stimdist','stimxdist','stimydist'});
disttab = disttab(ismember(disttab.stimcell,stimcells.cell),:);

% get column for distance from stim 

distpattab = outerjoin(pattab, disttab,'LeftKeys',{'cell'},'RightKeys',{'stimcell'}, 'Type','Left' ); % get all possible cell/stimcell pairs and their distances

trialstimcells = outerjoin(trialtab,pattab,'Keys',{'pattern','recid'}, 'Type','Left'); % get the cells stimmed in each trial
trialstimcells = trialstimcells(:,{'trial','cell','recid_trialtab'});
% tmp = join(stattab,trialstimcells,'LeftKeys',{'trial','recid'},'RightKeys',{'trial','recid_trialtab'});
% tmp.Properties.VariableNames{'cell_stattab'} = 'cell';
% tmp.Properties.VariableNames{ 'cell_trialstimcells'} =  'stimcell';
% tmp = tmp(tmp.stimcell>0,:);
% tmp = join(tmp,disttab,'LeftKeys',{'cell','stimcell'},'RightKeys',{'viewcell','stimcell'});
% diststats = grpstats(tmp, {'cell','recid'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});

tmp2 = outerjoin(disttab,pattab, 'LeftKeys',{'stimcell'},'RightKeys',{'cell'},'RightVariables',{'pattern','recid'}, 'Type', 'right'); % get the patterns containing any given stim cell
min_by_pattern = grpstats(tmp2, {'viewcell','pattern','recid'},{'min'},'DataVars','stimdist');
mindisttab = join(min_by_pattern,tmp2,'LeftKeys',{'viewcell','recid','pattern','min_stimdist'},'RightKeys',{'viewcell','recid','pattern','stimdist'});
stattab = join(stattab,mindisttab,'LeftKeys',{'cell','recid','pattern'},'RightKeys',{'viewcell','recid','pattern'},'RightVariables',{'min_stimdist','stimxdist','stimydist','stimcell'});
stattab.Properties.VariableNames{'min_stimdist'} = 'stimdist';
clear tmp tmp2 mindisttab trialstimcells;



% find bad cells based on min predff threshold
min_thresh = 20;
tmp = grpstats(stattab,"cell",{'mean','min','max'},'DataVars',{'dff','spontdff','premean'});
badcells = tmp.cell(tmp.min_premean < min_thresh);
clear tmp;

%% calculate stats for stim cells

cellstim_groups = findgroups(stattab.cell, stattab.cellstim);
cellstim_resp_groups = findgroups(stattab.cell, stattab.cellstim, stattab.resptype);

cellstats = grpstats(stattab, {'cell','recid'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});
cellstats_onestim = grpstats(stattab, {'cell', 'onestim'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});
stimstats = grpstats(stattab, {'cell', 'cellstim'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});
stimstats_resptype = grpstats(stattab, {'cell', 'cellstim', 'resptype'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});
stimstats_celltype = grpstats(stattab, { 'cellstim', 'celltype'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});
stimstats_celltype_resptype = grpstats(stattab, { 'cellstim', 'celltype', 'resptype'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});

stimstats_onestim = grpstats(stattab, {'cell', 'cellstim','onestim'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});
stimstats_onestim_resptype = grpstats(stattab, {'cell', 'cellstim', 'resptype','onestim'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});

% get correlation within groups
cellstats.dffstimcorr = splitapply(@(x,y) corr(x,y), double(stattab.respsnr), double(stattab.cellstim), findgroups(stattab.cell, stattab.recid));
cellstats_onestim.dffstimcorr = splitapply(@(x,y) corr(x,y), double(stattab.respsnr), double(stattab.cellstim), findgroups(stattab.cell, stattab.onestim));

% ks test statistic for cell-stimcell pairs - compare stimcell-stimmed distribution
% to distribution of all other trials minus those where the cell itself is
% stimmed

stimstats_paired = grpstats(stattab, { 'cell', 'stimcell'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});
ks_pval = zeros(size(stimstats_paired,1),1);
for irow=1:size(stimstats_paired,1)
    icell = stimstats_paired{irow,'cell'};
    istimcell = stimstats_paired{irow, 'stimcell'};
    %[~,ks_pval(irow)] = kstest2(stattab.dff(stattab.stimcell==istimcell & stattab.cell == icell), stattab.dff(stattab.stimcell ~= istimcell & stattab.stimcell ~= icell & stattab.cell == icell));
    [~,ks_pval(irow)] = kstest2(stattab.dff(stattab.stimcell==istimcell & stattab.cell == icell), stattab.spontdff(stattab.cell == icell));
end

stimstats_paired{:,'ks_pval'} = ks_pval;


%% calculate mean response curves for each cell


norm_by_f0 = true; % if true, responses are DF/F0 instead of DF

% set original movie or residuals here
%movmap = Fcent_list_bgsub; % for background -subtracted
movmap = Fcent_list; % for raw traces


% all the following are indexed by [index of cell in stimcells.cell, index of cell being measured ]
mean_resp = cell(nstimcells,ncells);
std_resp = cell(nstimcells,ncells);
mean_resp_nostim = cell(1,ncells);
std_resp_nostim = cell(1,ncells);



prestimindsfull = 1:(fullWindowPreSize - Omitpre);
poststimindsfull = poststiminds(1):size(movmap(1),1);

%mean_recids = onestim_recids; % movies to use
mean_recids = filetab.recid;

for istimcell=1:nstimcells
    stimcellnum = stimcells.cell(istimcell);
    tmpstat = stattab(:,{'trial','recid','cell','cellstim'});
    for recind=1:numel(mean_recids)
        recid = mean_recids(recind);
        Fc = movmap(recid);
        tmptrialtab = trialtab(trialtab.recid == recid,:);
        %inittrial = min(tmptrialtab.trial);
        inittrial = 1;
        stimtrialnums = tmpstat.trial(tmpstat.recid == recid &  tmpstat.cellstim & tmpstat.cell==stimcellnum);
        nostimtrialnums = tmpstat.trial(tmpstat.recid == recid  & ~tmpstat.cellstim & tmpstat.cell==stimcellnum);
        
        tmptrace = Fc(:,stimtrialnums - inittrial + 1, :); % first open index is frame in window, last open index is cell number
        tmptrace_nostim = Fc(:,nostimtrialnums - inittrial + 1, :); % first open index is frame in window, last open index is cell number
        
        if norm_by_f0
           tmptrace = ( tmptrace([prestimindsfull,poststimindsfull],:,:) - mean(tmptrace(prestimindsfull,:,:),1) ) ./ mean(tmptrace(prestimindsfull,:,:),1) ;   
           tmptrace_nostim = ( tmptrace_nostim([prestimindsfull,poststimindsfull],:,:) - mean(tmptrace_nostim(prestimindsfull,:,:),1) ) ./ mean(tmptrace_nostim(prestimindsfull,:,:),1) ;   
        else
            tmptrace = ( tmptrace([prestimindsfull,poststimindsfull],:,:) - mean(tmptrace(prestimindsfull,:,:),1) )  ;   
            tmptrace_nostim = ( tmptrace_nostim([prestimindsfull,poststimindsfull],:,:) - mean(tmptrace_nostim(prestimindsfull,:,:),1) ) ;
        end
        
        if recind==1
            traces = tmptrace;
            traces_nostim = tmptrace_nostim;
        else
            traces = cat(2,traces, tmptrace);
            traces_nostim = cat(2,traces_nostim, tmptrace_nostim);
        end
    end
    for meascellnum=1:ncells
        mean_resp{istimcell,meascellnum} = mean(traces(:,:,meascellnum),2);
        std_resp{istimcell,meascellnum} = std(traces(:,:,meascellnum),0,2) / sqrt(size(traces,2));
        mean_resp_nostim{1,meascellnum} = mean(traces_nostim(:,:,meascellnum),2);
        std_resp_nostim{1,meascellnum} = std(traces_nostim(:,:,meascellnum),0,2) / sqrt(size(traces_nostim,2));
    end
end

clear tmptrace tmptrace_nostim traces_nostim
% %% plots
% 
% % plot stim cell responses
% selfstiminds = sub2ind(size(mean_resp),1:nstimcells,stimcells.cell');
% 
% 
% tmp = cat(2,mean_resp{selfstiminds});
% figure; hold on; 
% myaxes = [-inf inf -.05 .1];
% plot( tmp); axis(myaxes); h=plot([numel(prestimindsfull) numel(prestimindsfull)], myaxes([3 4]),'LineStyle','-','Color','r', 'LineWidth',2);
% hold off;
% 
% 
% 
% tmp = cat(2,mean_resp{stimcells.cell});
% figure; hold on; 
% myaxes = [-inf inf -.1 .15];
% plot( tmp); axis(myaxes); h=plot([numel(prestimindsfull) numel(prestimindsfull)], myaxes([3 4]),'LineStyle','-','Color','r','LineWidth',2); 
% hold off;

%% plot cell stim profiles (targeted / nontargeted)

%istimcell = 1;
%cellnum = stimcells.cell(istimcell);
iresp = 4;
cellnum = responders(iresp);
istimcell = find(stimcells.cell==cellnum);

err_include = true;

 figure;
 gaps = [.1 .04];
 marg_w = .1;
 marg_h = .1;
 
axislim =  [-inf inf -.03 .085];
 
graycolor = [.5 .5 .5];

% histogram for red cells during responder stim
ax = subtightplot(1,2,1,gaps,marg_h, marg_w);
mean_vals = mean_resp{istimcell,cellnum}; sem_vals = std_resp{istimcell,cellnum};
%mean_vals = mean_resp{14,11}; sem_vals = std_resp{14,11};
plotwindow = [prestimindsfull poststimindsfull] - fullWindowPreSize;
li = plot(plotwindow, mean_vals , 'Color','r');
li.Color = graycolor;
if err_include
    hold on;
    er = shadedErrorBar(plotwindow, mean_vals,sem_vals);
    er.LineStyle = 'none';
    er.Color = graycolor;
    hold off;
end
ylabel('DF/F_0');
axis(axislim);
title(sprintf('targeted stim response: cell %d', cellnum));
xlabel('Frames relative to stim time');

% histogram for red cells during responder stim
ax = subtightplot(1,2,2,gaps,marg_h, marg_w);
nonstiminds = setdiff(1:size(mean_resp,1),istimcell);
mean_vals = mean_resp_nostim{1,cellnum}; sem_vals = std_resp_nostim{1,cellnum};
plotwindow = [prestimindsfull poststimindsfull] - fullWindowPreSize;
li = plot(plotwindow, mean_vals , 'Color', 'r');
li.Color = graycolor;
if err_include
    hold on;
    er = shadedErrorBar(plotwindow, mean_vals,sem_vals);
    er.LineStyle = 'none';
    er.Color = graycolor;
    hold off;
end
axis(axislim);
set(gca,'YTickLabel',[]);
title('non-targeted stim response');

%% stim responses (nonstim cells)

cellnum = 1;
err_include = 1;

 figure;
 gaps = [.1 .04];
 marg_w = .1;
 marg_h = .1;
 
axislim =  [-inf inf -.05 .085];
 
graycolor = [.5 .5 .5];

mean_vals = mean_resp_nostim{1,cellnum}; sem_vals = std_resp_nostim{1,cellnum};
plotwindow = [prestimindsfull poststimindsfull] - fullWindowPreSize;
li = plot(plotwindow, mean_vals , 'Color','r');
li.Color = graycolor;
if err_include
    hold on;
    er = shadedErrorBar(plotwindow, mean_vals,sem_vals);
    er.LineStyle = 'none';
    er.Color = graycolor;
    hold off;
end
ylabel('DF/F_0');
axis(axislim);
title(sprintf('non-targeted stim response: cell %d', cellnum));
xlabel('Frames relative to stim time');


%% cross correlation matrix between stimulated cell indicator and observed cell dff

[grp, id] = findgroups(stattab.cell);
corrmat = zeros(nstimcells,ncells);

for istimcell=1:numel(stimcells.cell)
    stimcellnum = stimcells.cell(istimcell);
    stimtrials = unique(stattab.trial(stattab.cell == stimcellnum & stattab.cellstim));
    corrmat(istimcell,id) = splitapply(@corr ,stattab.dff, ismember(stattab.trial,stimtrials) , grp); 
end

corrmatref = zeros(size(corrmat)); corrmatref(sub2ind(size(corrmatref),1:size(corrmatref,1),stimcells.cell')) = 1; % indicated how stimcells and cells match
figure; imagesc(corrmat); colorbar;caxis([-.08 .08]); colormap(mymap);
figure; imagesc(corrmatref); caxis([-1 1]);





%% dff vs distance


plotdiff = true;
adjscale=false;


% dff

 figure;
 gaps = [.1 .04];
 marg_w = .1;
 marg_h = .1;
 

bin_partition = 0:20:400;
myaxes = [-20 bin_partition(end) -.01 .05];
 
graycolor = [.5 .5 .5];

% histogram for red cells during responder stim
indicator = stattab.celltype==1 & ismember(stattab.stimcell,find(isresponder)) &...
    stattab.stimdist <= max(bin_partition) & ... 
    stattab.stimdist >= min(bin_partition) & ~ismember(stattab.cell,badcells) & ...
     true ;%(ismember(stattab.cell,responders)  ) ;
dists = stattab.stimdist(indicator);
dffs = stattab.dff(indicator);
spontdffs = stattab.randdff(indicator);
[binnums, edges] = discretize(dists,bin_partition);
[binmeans,binstds, bincounts] = grpstats(dffs,binnums,{'mean','std','numel'});
usedbins = sort(unique(binnums));
binstds = binstds ./ sqrt(bincounts);
[binmeansspont,binstdsspont,~] = grpstats(spontdffs,binnums,{'mean','std','numel'});
binstdsspont = binstdsspont ./ sqrt(bincounts);
ploterrbound =  binstdsspont;
empty_edges = setdiff(1:numel(edges),usedbins);
edges(empty_edges) = [];
hold on;
li = plot(edges,binmeans, 'o', 'MarkerEdgeColor', 'red', 'MarkerFaceColor','red', 'MarkerSize',3);
li.Color = 'red';
er = errorbar(edges, binmeans,binstds, binstds);
er.Color = 'red';
adjmax = max(binmeans);

ylabel('DF/F_0');
title('stim response');
axis(myaxes);
xlabel('Distance from nearest stim site (pixels)');


% histogram for non-red cells
indicator = ismember(stattab.stimcell,find(isresponder)) & stattab.celltype==0 & stattab.stimdist <= max(bin_partition) & stattab.stimdist >= min(bin_partition) & ~ismember(stattab.cell,badcells)  ;
dists = stattab.stimdist(indicator);
dffs = stattab.dff(indicator);
spontdffs = stattab.randdff(~ismember(stattab.cell,badcells));
[binnums, edges] = discretize(dists,bin_partition);
[spontbinnums, spontedges] = discretize(stattab.stimdist(~ismember(stattab.cell,badcells)),bin_partition);
[binmeans,binstds, bincounts] = grpstats(dffs,binnums,{'mean','std','numel'});
usedbins = sort(unique(binnums));
binstds = binstds ./ sqrt(bincounts);
[binmeansspont,binstdsspont,bincountsspont] = grpstats(spontdffs,spontbinnums,{'mean','std','numel'});
binstdsspont = binstdsspont ./ sqrt(bincountsspont);
ploterrbound = binstdsspont;
empty_edges = setdiff(1:numel(edges),usedbins);
edges(empty_edges) = [];
empty_spontedges = setdiff(1:numel(spontedges),spontbinnums);
spontedges(empty_spontedges) = [];
hold on;
if plotdiff
plot(spontedges, binmeansspont, 'x', 'MarkerEdgeColor', 'black', 'MarkerFaceColor','black' );
er = errorbar(spontedges, binmeansspont,ploterrbound, ploterrbound);
er.Color = 'black';
end
li = plot(edges,binmeans, 'o', 'MarkerEdgeColor', 'green', 'MarkerFaceColor','green', 'MarkerSize',3);
li.Color = 'green';
er = errorbar(edges, binmeans,binstds, binstds);
er.Color = 'green';
grid on;
hold off;
axis square;

if adjscale
    yticks([0 .25 .5 .75 1]*adjmax);
    yticklabels([0 .25 .5 .75 1]);
end


%% activation / inhibition probability vs distance, based on dff

 
 gaps = [.1 .08];
 marg_w = .1;
 marg_h = .09;
 
 
bin_partition = -30:40:400;
myaxes = [-20 bin_partition(end) -.2 .6];

activate_thresh = .03; % dff threshold for activation
inhibit_thresh = -.03;

activate_data = stattab.dff > activate_thresh;
inhibit_data = stattab.dff < inhibit_thresh;
spont_activate_data = stattab.spontdff > activate_thresh;
spont_inhibit_data = stattab.spontdff < inhibit_thresh;


figure;
% activation for red cells during responder stim
ax = subtightplot(2,2,1,gaps,marg_h, marg_w);
indicator = stattab.celltype==1 & ismember(stattab.stimcell,find(isresponder)) &...
    stattab.stimdist <= max(bin_partition) & ... 
    stattab.stimdist >= min(bin_partition) & ~ismember(stattab.cell,badcells) & ...
     true ;%(ismember(stattab.cell,responders)  ) ;
dists = stattab.stimdist(indicator);
dffs = activate_data(indicator);
spontdffs = spont_activate_data(indicator);
[binnums, edges] = discretize(dists,bin_partition);
[binmeans,binstds, bincounts] = grpstats(dffs,binnums,{'mean','std','numel'});
usedbins = sort(unique(binnums));
binstds = binstds ./ sqrt(bincounts);
[binmeansspont,binstdsspont,~] = grpstats(spontdffs,binnums,{'mean','std','numel'});
binstdsspont = binstdsspont ./ sqrt(bincounts);
ploterrbound = 2 * ( binstdsspont + binstds ) / sqrt(2);
empty_edges = setdiff(1:numel(edges),usedbins);
edges(empty_edges) = [];
hold on; edges = edges + 10;
bar(edges, binmeans - binmeansspont);
li = plot(edges,binmeans);
li.Color = [.5 .5 .5];
er = errorbar(edges, binmeans - binmeansspont,ploterrbound, ploterrbound);
er.LineStyle = 'none';
er.Color = 'red';
hold off;
ylabel('activation probability');
title('red label stim response');
axis(myaxes);

% activation for non-red cells during responder stim
ax = subtightplot(2,2,2,gaps,marg_h, marg_w);
indicator = stattab.celltype==0 & ismember(stattab.stimcell,find(isresponder)) &...
    stattab.stimdist <= max(bin_partition) & ... 
    stattab.stimdist >= min(bin_partition) & ~ismember(stattab.cell,badcells) & ...
     true ;%(ismember(stattab.cell,responders)  ) ;
dists = stattab.stimdist(indicator);
dffs = activate_data(indicator);
spontdffs = spont_activate_data(indicator);
[binnums, edges] = discretize(dists,bin_partition);
[binmeans,binstds, bincounts] = grpstats(dffs,binnums,{'mean','std','numel'});
binstds = binstds ./ sqrt(bincounts);
usedbins = sort(unique(binnums));
[binmeansspont,binstdsspont,~] = grpstats(spontdffs,binnums,{'mean','std','numel'});
binstdsspont = binstdsspont ./ sqrt(bincounts);
ploterrbound = 2 * ( binstdsspont + binstds ) / sqrt(2);
empty_edges = setdiff(1:numel(edges),usedbins);
edges(empty_edges) = [];
hold on; edges = edges + 10;
bar(edges, binmeans - binmeansspont);
li = plot(edges,binmeans);
li.Color = [.5 .5 .5];
er = errorbar(edges, binmeans - binmeansspont,ploterrbound, ploterrbound);
er.LineStyle = 'none';
er.Color = 'red';
hold off;
title('non-red label stim response');
axis(myaxes);

% inhibition for red cells during responder stim
ax = subtightplot(2,2,3,gaps,marg_h, marg_w);
indicator = ismember(stattab.stimcell,find(isresponder)) & stattab.celltype==1 & stattab.stimdist <= max(bin_partition) & ... 
    stattab.stimdist >= min(bin_partition) & ~ismember(stattab.cell,badcells) & ...
     true; %(ismember(stattab.cell,responders) ) ;
dists = stattab.stimdist(indicator);
dffs = inhibit_data(indicator);
spontdffs = spont_inhibit_data(indicator);
[binnums, edges] = discretize(dists,bin_partition);
[binmeans,binstds, bincounts] = grpstats(dffs,binnums,{'mean','std','numel'});
binstds = binstds ./ sqrt(bincounts);
[binmeansspont,binstdsspont,~] = grpstats(spontdffs,binnums,{'mean','std','numel'});
binstdsspont = binstdsspont ./ sqrt(bincounts);
ploterrbound = 2 * ( binstdsspont + binstds ) / sqrt(2);
empty_edges = setdiff(1:numel(edges),binnums);
edges(empty_edges) = [];
hold on; edges = edges + 10;
bar(edges, binmeans - binmeansspont);
li = plot(edges,binmeans);
li.Color = [.5 .5 .5];
er = errorbar(edges, binmeans - binmeansspont,ploterrbound, ploterrbound);
er.LineStyle = 'none';
er.Color = 'red';
hold off;
ylabel('inhibition probability');
xlabel('distance from nearest stim site (pixels)');
axis(myaxes);


% inhibition for non-red cells during responder stim
ax = subtightplot(2,2,4,gaps,marg_h, marg_w);
indicator = ismember(stattab.stimcell,find(isresponder)) & stattab.celltype==0 & stattab.stimdist <= max(bin_partition) & stattab.stimdist >= min(bin_partition) & ~ismember(stattab.cell,badcells);;
dists = stattab.stimdist(indicator);
dffs = inhibit_data(indicator);
spontdffs = spont_inhibit_data(indicator);
[binnums, edges] = discretize(dists,bin_partition);
[binmeans,binstds, bincounts] = grpstats(dffs,binnums,{'mean','std','numel'});
binstds = binstds ./ sqrt(bincounts);
[binmeansspont,binstdsspont,~] = grpstats(spontdffs,binnums,{'mean','std','numel'});
binstdsspont = binstdsspont ./ sqrt(bincounts);
ploterrbound = 2 * ( binstdsspont + binstds ) / sqrt(2);
empty_edges = setdiff(1:numel(edges),binnums);
edges(empty_edges) = [];
hold on; edges = edges + 10;
bar(edges, binmeans - binmeansspont);
li = plot(edges,binmeans);
li.Color = [.5 .5 .5];
er = errorbar(edges, binmeans - binmeansspont,ploterrbound, ploterrbound);
er.LineStyle = 'none';
er.Color = 'red';
hold off;
axis(myaxes);





%% get spatial activation map (red cells)
blursigma = 5; % blurring filter
mapsize = 2*[height width];
mapcenter = round(mapsize/2);
finalzoom = 1.7;




clims_prob = [-.3,.3];
clims = [-.1 .1];
circle_cells = false; %if true, use approximating circles for cell shapes

% check the recording type here...
spatial_tab = grpstats(stattab(ismember(stattab.stimcell,find(isresponder)) & stattab.celltype==1 & ~ismember(stattab.cell,badcells) ,:),{'cell','stimcell','resptype'},{'mean'},'DataVars',{'dff','stimxdist','stimydist'});
spatial_spots = sub2ind(mapsize, mapcenter(1)-round(spatial_tab.mean_stimydist), mapcenter(2)-round(spatial_tab.mean_stimxdist));

    


spatial_map_excite = zeros(mapsize);
spatial_map_inhibit = zeros(mapsize);
spatial_map_dff = zeros(mapsize);
spatial_map_base = zeros(mapsize); %stores "denominator image" for map (local number of cells)


for i=1:size(spatial_tab,1)
    icell = spatial_tab{i,'cell'};
    if circle_cells
        update_spot = spatial_spots(i);
    else
        [y,x] = ind2sub([height width], maskinds{icell});
        y = round(y - yavg(icell) - spatial_tab.mean_stimydist(i) + mapcenter(1));
        x = round(x - xavg(icell) - spatial_tab.mean_stimxdist(i) + mapcenter(2));
        update_spot = sub2ind(mapsize,y,x);
    end
    if spatial_tab.resptype(i)==1
        spatial_map_excite(update_spot) =  spatial_map_excite(update_spot) + spatial_tab.GroupCount(i);
    elseif spatial_tab.resptype(i)==-1
       spatial_map_inhibit(update_spot) =  spatial_map_inhibit(update_spot) + spatial_tab.GroupCount(i);
    end
    spatial_map_dff(update_spot) =  spatial_map_dff(update_spot) + spatial_tab.mean_dff(i)*spatial_tab.GroupCount(i);
    spatial_map_base(update_spot) = spatial_map_base(update_spot) + spatial_tab.GroupCount(i);
end

    
spatial_map_excite =  imgaussfilt(spatial_map_excite, blursigma,'FilterSize',4*ceil(2*blursigma)+1.);
spatial_map_inhibit = imgaussfilt(spatial_map_inhibit, blursigma,'FilterSize',4*ceil(2*blursigma)+1.);
spatial_map_dff = imgaussfilt(spatial_map_dff, blursigma,'FilterSize',4*ceil(2*blursigma)+1.);
spatial_map_base = imgaussfilt(spatial_map_base, blursigma,'FilterSize',4*ceil(2*blursigma)+1.);

spatial_map_base(spatial_map_base == 0) = inf;
spatial_map_excite = spatial_map_excite ./ spatial_map_base;
spatial_map_inhibit = spatial_map_inhibit ./ spatial_map_base;
spatial_map_dff = spatial_map_dff ./ spatial_map_base ;


spatial_map_difference = spatial_map_excite - spatial_map_inhibit;


% center circle for plotting
th = 0:.1:2*pi+.1;
r = 20;
xunit = r * cos(th) + mapcenter(2);
yunit = r * sin(th) + mapcenter(1);

figure;
subtightplot(2,4,1);
hold on;
imagesc(spatial_map_excite, clims_prob); colorbar; colormap(mymap);
plot(xunit, yunit); title('red cell spatial excitation map');
xlim([mapcenter(2)-round(width/(2*finalzoom)),mapcenter(2)+round(width/(2*finalzoom))]); ylim([mapcenter(1)-round(height/(2*finalzoom)),mapcenter(1)+round(height/(2*finalzoom))]);
hold off;
subtightplot(2,4,2);
hold on;
imagesc(-spatial_map_inhibit, clims_prob); colorbar;
plot(xunit, yunit); title('red cell spatial inhibition map');
xlim([mapcenter(2)-round(width/(2*finalzoom)),mapcenter(2)+round(width/(2*finalzoom))]); ylim([mapcenter(1)-round(height/(2*finalzoom)),mapcenter(1)+round(height/(2*finalzoom))]);
set(gca,'YTickLabel',[]);
hold off;
subtightplot(2,4,3);
hold on;
imagesc(spatial_map_difference, clims_prob); colorbar;
plot(xunit, yunit); title('red cell spatial difference prob. map');
xlim([mapcenter(2)-round(width/(2*finalzoom)),mapcenter(2)+round(width/(2*finalzoom))]); ylim([mapcenter(1)-round(height/(2*finalzoom)),mapcenter(1)+round(height/(2*finalzoom))]);
set(gca,'YTickLabel',[]);
hold off;
subtightplot(2,4,4);
hold on;
imagesc(spatial_map_dff, clims); h=colorbar;colormap(mymap);
ylabel(h,'local mean DF/F0');
plot(xunit, yunit); title('red cell spatial dff map');
xlim([mapcenter(2)-round(width/(2*finalzoom)),mapcenter(2)+round(width/(2*finalzoom))]); ylim([mapcenter(1)-round(height/(2*finalzoom)),mapcenter(1)+round(height/(2*finalzoom))]);
set(gca,'YTickLabel',[]);
hold off;
set(gca,'Layer','top');

% axis square
% xticks(512+(-200:50:200))
% xticklabels((-200:50:200))
% yticks(512+(-200:50:200))
% yticklabels((-200:50:200))
% xlabel('Displacement from stim site (pixels)');

%% get spatial activation map (nonred cells)
blursigma = 5; % blurring filter
mapsize = 2*[height width];
mapcenter = round(mapsize/2);
finalzoom = 2;


nmap = 128;
cmin = .2;
mymap = vertcat([linspace(cmin,1,nmap)' linspace(cmin,1,nmap)' ones(nmap,1)], [ones(nmap,1) linspace(1,cmin,nmap)' linspace(1,cmin,nmap)']); % blue->white->red


clims_prob = [-.3,.3];
clims = [-.05 .05];
circle_cells = false; %if true, use approximating circles for cell shapes

% check the recording type here...
spatial_tab = grpstats(stattab(ismember(stattab.stimcell,find(isresponder)) & stattab.celltype==0 & ~ismember(stattab.cell,badcells) ,:),{'cell','stimcell','resptype'},{'mean'},'DataVars',{'dff','stimxdist','stimydist'});
spatial_spots = sub2ind(mapsize, mapcenter(1)-round(spatial_tab.mean_stimydist), mapcenter(2)-round(spatial_tab.mean_stimxdist));


spatial_map_excite = zeros(mapsize);
spatial_map_inhibit = zeros(mapsize);
spatial_map_dff = zeros(mapsize);
spatial_map_base = zeros(mapsize); %stores "denominator image" for map (local number of cells)


for i=1:size(spatial_tab,1)
    icell = spatial_tab{i,'cell'};
    if circle_cells
        update_spot = spatial_spots(i);
    else
        [y,x] = ind2sub([height width], maskinds{icell});
        y = round(y - yavg(icell) - spatial_tab.mean_stimydist(i) + mapcenter(1));
        x = round(x - xavg(icell) - spatial_tab.mean_stimxdist(i) + mapcenter(2));
        update_spot = sub2ind(mapsize,y,x);
    end
    if spatial_tab.resptype(i)==1
        spatial_map_excite(update_spot) =  spatial_map_excite(update_spot) + spatial_tab.GroupCount(i);
    elseif spatial_tab.resptype(i)==-1
       spatial_map_inhibit(update_spot) =  spatial_map_inhibit(update_spot) + spatial_tab.GroupCount(i);
    end
    spatial_map_dff(update_spot) =  spatial_map_dff(update_spot) + spatial_tab.mean_dff(i)*spatial_tab.GroupCount(i);
    spatial_map_base(update_spot) = spatial_map_base(update_spot) + spatial_tab.GroupCount(i);
end

    
spatial_map_excite =  imgaussfilt(spatial_map_excite, blursigma,'FilterSize',4*ceil(2*blursigma)+1.);
spatial_map_inhibit = imgaussfilt(spatial_map_inhibit, blursigma,'FilterSize',4*ceil(2*blursigma)+1.);
spatial_map_dff = imgaussfilt(spatial_map_dff, blursigma,'FilterSize',4*ceil(2*blursigma)+1.);
spatial_map_base = imgaussfilt(spatial_map_base, blursigma,'FilterSize',4*ceil(2*blursigma)+1.);

spatial_map_base(spatial_map_base == 0) = inf;
spatial_map_excite = spatial_map_excite ./ spatial_map_base;
spatial_map_inhibit = spatial_map_inhibit ./ spatial_map_base;
spatial_map_dff = spatial_map_dff ./ spatial_map_base ;


spatial_map_difference = spatial_map_excite - spatial_map_inhibit;


% center circle for plotting
th = 0:.1:2*pi+.1;
r = 20;
xunit = r * cos(th) + mapcenter(2);
yunit = r * sin(th) + mapcenter(1);

figure;
subtightplot(2,4,1);
hold on;
imagesc(spatial_map_excite, clims_prob); colorbar; colormap(mymap);
plot(xunit, yunit); title('nonred cell spatial excitation map');
xlim([mapcenter(2)-round(width/(2*finalzoom)),mapcenter(2)+round(width/(2*finalzoom))]); ylim([mapcenter(1)-round(height/(2*finalzoom)),mapcenter(1)+round(height/(2*finalzoom))]);
hold off;
subtightplot(2,4,2);
hold on;
imagesc(-spatial_map_inhibit, clims_prob); colorbar;
plot(xunit, yunit); title('nonred cell spatial inhibition map');
xlim([mapcenter(2)-round(width/(2*finalzoom)),mapcenter(2)+round(width/(2*finalzoom))]); ylim([mapcenter(1)-round(height/(2*finalzoom)),mapcenter(1)+round(height/(2*finalzoom))]);
set(gca,'YTickLabel',[]);
hold off;
subtightplot(2,4,3);
hold on;
imagesc(spatial_map_difference, clims_prob); colorbar;
plot(xunit,yunit); title('nonred cell spatial difference prob. map');
xlim([mapcenter(2)-round(width/(2*finalzoom)),mapcenter(2)+round(width/(2*finalzoom))]); ylim([mapcenter(1)-round(height/(2*finalzoom)),mapcenter(1)+round(height/(2*finalzoom))]);
set(gca,'YTickLabel',[]);
hold off;
subtightplot(2,4,4);
hold on;
imagesc(spatial_map_dff, clims); colorbar;
plot(xunit,yunit); title('nonred cell spatial dff map');
xlim([mapcenter(2)-round(width/(2*finalzoom)),mapcenter(2)+round(width/(2*finalzoom))]); ylim([mapcenter(1)-round(height/(2*finalzoom)),mapcenter(1)+round(height/(2*finalzoom))]);
set(gca,'YTickLabel',[]);
hold off;


%% full field activation map (individual stim)

blursigma = .1; % blurring filter
mapsize = [height width];
ks_resp_only = false; % only show cells with appropriate ks response threshold
ks_resp_thresh = .05;

load(multi_files{2},'stattab','isresponder','badcells', 'responders', 'maskinds','stimcells', 'stimstats_paired','xavg', 'yavg', 'maskinds');



plotcells = responders; 
export_dir = ['C:\Users\bnste\Downloads\JG16271\JG16271_stim_response_plots_10']; % leave empty to not save

clims_prob = [-.3,.3];
clims = [-.1 .1];
circle_cells = false; %if true, use approximating circles for cell shapes


cell_type = 2; % 0 for nonred, 1 for red, 2 for all


for icellnum=1:numel(plotcells)


cellnum = plotcells(icellnum);
stimind = find(stimcells.cell,cellnum);
 

% check the recording type here...
indicator = stattab.stimcell == cellnum & ismember(stattab.stimcell,find(isresponder))  & ismember(stattab.recid, onestim_recids) & stattab.celltype >= 0 & ~ismember(stattab.cell,badcells);
if cell_type < 2
    indicator = indicator & stattab.celltype==cell_type; 
end
spatial_tab = grpstats(stattab( indicator,:),{'cell','stimcell','resptype'},{'mean'},'DataVars',{'dff','stimxdist','stimydist'});

if ks_resp_only
   spatial_tab = spatial_tab(ismember(spatial_tab.cell, stimstats_paired.cell(stimstats_paired.ks_pval < ks_resp_thresh & stimstats_paired.stimcell == cellnum))  ,:) ;
end
    


spatial_map_excite = zeros(mapsize);
spatial_map_inhibit = zeros(mapsize);
spatial_map_dff = zeros(mapsize);
spatial_map_base = zeros(mapsize); %stores "denominator image" for map (local number of cells)


for i=1:size(spatial_tab,1)
    icell = spatial_tab{i,'cell'};
 
    [y,x] = ind2sub([height width], maskinds{icell});
    update_spot = sub2ind(mapsize,y,x);

    if spatial_tab.resptype(i)==1
        spatial_map_excite(update_spot) =  spatial_map_excite(update_spot) + spatial_tab.GroupCount(i);
    elseif spatial_tab.resptype(i)==-1
       spatial_map_inhibit(update_spot) =  spatial_map_inhibit(update_spot) + spatial_tab.GroupCount(i);
    end
    spatial_map_dff(update_spot) =  spatial_map_dff(update_spot) + spatial_tab.mean_dff(i)*spatial_tab.GroupCount(i);
    spatial_map_base(update_spot) = spatial_map_base(update_spot) + spatial_tab.GroupCount(i);
end

spatial_map_base(spatial_map_base == 0) = inf;

spatial_map_dff = imgaussfilt(spatial_map_dff ./ spatial_map_base, blursigma,'FilterSize',4*ceil(2*blursigma)+1.);


% stim circle for plotting
th = 0:.1:2*pi;
r = 5;
xunit = r * cos(th) + xavg(cellnum);
yunit = r * sin(th) + yavg(cellnum);




figure;
hold on;
imagesc(spatial_map_dff, clims); 
h=colorbar; set(get(h,'label'), 'string', 'mean DF/F stim response'); 
colormap(mymap)
xlim([0 mapsize(2)]); ylim([0 mapsize(1)]); 
set(gca,'YDir','reverse');
title(sprintf('spatial dff map: stim cell %d',cellnum));
xlabel('pixels');
axis square;

plot(xunit,yunit);

c = redblue_log;
colormap(c);


hold off;

if numel(export_dir)>0
    if ks_resp_only
        saveas(gcf,fullfile(export_dir,sprintf('globalrespmap_thresh_%05d.eps', cellnum)),'epsc');
    else
        saveas(gcf,fullfile(export_dir,sprintf('globalrespmap_%05d.eps', cellnum)),'epsc');
    end
end
end
%% sample responses to go with above full field maps


%load(multi_files{1},'stattab','isresponder','badcells', 'responders', 'maskinds','stimcells', 'stimstats_paired','xavg', 'yavg', 'maskinds','mean_resp');
fps=30;

stimind = 13;
cellnum = stimcells.cell(stimind);
viewcell = 52;
step = 100; % ms
labelfrac = 5;

err_include = true;

 figure;

axislim =  [-inf inf -.04 .085];
 
graycolor = [.5 .5 .5];

mean_vals = mean_resp{stimind,viewcell}; sem_vals = std_resp{stimind,viewcell};
%mean_vals = mean_resp{14,11}; sem_vals = std_resp{14,11};
plotwindow = [prestimindsfull poststimindsfull] - fullWindowPreSize;
plotms = 1000 * ((1:numel(plotwindow)) - fullWindowPreSize+2) / fps;
gapframesize = min(poststimindsfull) - max(prestimindsfull) - 1;

li = plot(plotms,mean_vals , 'Color','r');
li.Color = graycolor;
if err_include
    hold on;
    er = shadedErrorBar(plotms, mean_vals,sem_vals);
    er.LineStyle = 'none';
    er.Color = graycolor;
    hold off;
end
ylabel('DF/F_0');


xtickangle(0);
xline(0,'--r','Linewidth',1);
xticks( [(min(xticks):step:-1) (0:step:max(xticks))]);
xtickformat

axis(axislim);
title(sprintf('cell %d response to cell %d stim', viewcell, cellnum));
xlabel('Time relative to stim period (ms)');


% lab=xticklabels;
% inds = find( mod(round(xticks/step),labelfrac)~=labelfrac-1);
% lab(inds) = {''};
% xticklabels(lab);


%% combining dff vs dist plots between mice

plotspont = true;
adjscale = false;

clear binstats*

% dff
 
 figure;
 gaps = [.1 .04];
 marg_w = .1;
 marg_h = .1;
 

bin_partition = 0:20:400;
myaxes = [-20 bin_partition(end) -.01 .05];
 
graycolor = [.5 .5 .5];
fulldffs = cell(nmultifiles,1);
fulldists = cell(nmultifiles,1);
fullspontdffs = cell(nmultifiles,1);
fullspontdists = cell(nmultifiles,1);
fullspontdists2 = cell(nmultifiles,1);
grpstats_red = cell(nmultifiles,1);
grpstats_nonred = cell(nmultifiles,1);
grpstats_all = cell(nmultifiles,1);


% histogram for red cells during responder stim

for iii=1:nmultifiles
    
    load(multi_files{iii},'stattab','isresponder','badcells', 'responders', 'maskinds');
     indicator = ismember(stattab.stimcell,find(isresponder)) & stattab.celltype==1 & stattab.stimdist <= max(bin_partition) & ... 
         stattab.stimdist >= min(bin_partition) & ~ismember(stattab.cell,badcells) ;%& (ismember(stattab.cell,responders)  ) ;
          
      
    fulldists{iii} = stattab.stimdist(indicator);
    fulldffs{iii} = stattab.dff(indicator);
    
    % get count statistics for red cells
    [binnums, ~] = discretize(stattab.stimdist(indicator),bin_partition);
    tmp = stattab(indicator,{'cell','stimcell'}); tmp.binnum = binnums;
    tmp = groupcounts(tmp, {'cell','stimcell','binnum'},'IncludeEmptyGroups', true );
    tmp.expid = repmat(iii, [size(tmp,1) 1]);
    grpstats_red{iii} = tmp;
    
     
    fullspontdffs{iii} = stattab.randdff(indicator);
    fullspontdists{iii} = stattab.stimdist(indicator);
    
end
dffs1 = cat(1,fulldffs{:}); dists = cat(1,fulldists{:}); spontdffs1 = cat(1,fullspontdffs{:}); spontdists1 = cat(1,fullspontdists{:});
[binnums1, edges] = discretize(dists,bin_partition);
[binmeans1,binstds1, bincounts1] = grpstats(dffs1,binnums1,{'mean','std','numel'});
binstds1 = binstds1 ./ sqrt(bincounts1);
empty_edges = setdiff(1:numel(edges),binnums1);
edges(empty_edges) = [];
[spontbinnums1, spontedges1] = discretize(spontdists1,bin_partition);
[binmeansspont1,binstdsspont1,bincountsspont1] = grpstats(spontdffs1,spontbinnums1,{'mean','std','numel'});
hold on;
li = plot(edges,binmeans1, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor','r', 'MarkerSize',3);
li.Color = 'r';
er = errorbar(edges, binmeans1,binstds1, binstds1);
er.Color = 'r';
hold off;
ylabel('DF/F_0');
title('red label stim response');
axis(myaxes);
xlabel('Distance from nearest stim site (microns)');
adjmax = max(binmeans1);


% non-red cells
for iii=1:nmultifiles
    load(multi_files{iii},'stattab','isresponder','badcells', 'responders', 'maskinds');
    indicator = ismember(stattab.stimcell,find(isresponder)) & stattab.celltype==0 & stattab.stimdist <= max(bin_partition) & ... 
        stattab.stimdist >= min(bin_partition) & ~ismember(stattab.cell,badcells)  ;
    fulldists{iii} = stattab.stimdist(indicator);
    fulldffs{iii} = stattab.dff(indicator);
    
    
    % get count statistics for nonred cells
    [binnums, ~] = discretize(stattab.stimdist(indicator),bin_partition);
    tmp = stattab(indicator,{'cell','stimcell'}); tmp.binnum = binnums;
    tmp = groupcounts(tmp, {'cell','stimcell','binnum'},'IncludeEmptyGroups', true );
    tmp.expid = repmat(iii, [size(tmp,1) 1]);
    grpstats_nonred{iii} = tmp;
    
    
    
    
    indicator = ismember(stattab.stimcell,find(isresponder)) & stattab.stimdist <= max(bin_partition) & ... 
    stattab.stimdist >= min(bin_partition) & ~ismember(stattab.cell,badcells) & stattab.celltype >= 0;% & ( ismember(stattab.cell,responders) | stattab.celltype == 0  ); % only include responder red cells
    indicator2 = ismember(stattab.stimcell,find(isresponder)) & stattab.stimdist <= max(bin_partition) & ... 
    stattab.stimdist >= min(bin_partition) & ~ismember(stattab.cell,badcells) & stattab.celltype==0  ;
    fullspontdists2{iii} = stattab.stimdist(indicator2);
    fullspontdffs2{iii} = stattab.randdff(indicator2);
    fullspontdists{iii} = stattab.stimdist(indicator);
    fullspontdffs{iii} = stattab.randdff(indicator);
    
    
    alldists{iii} = stattab.stimdist(indicator);
    alldffs{iii} = stattab.dff(indicator);
    
     % get count statistics for all cells
    [binnums, ~] = discretize(stattab.stimdist(indicator),bin_partition);
    tmp = stattab(indicator,{'cell','stimcell'}); tmp.binnum = binnums;
    tmp = groupcounts(tmp, {'cell','stimcell','binnum'},'IncludeEmptyGroups', true );
    tmp.expid = repmat(iii, [size(tmp,1) 1]);
    grpstats_all{iii} = tmp;
    
end
dffs2 = cat(1,fulldffs{:}); dists = cat(1,fulldists{:});spontdffs2 = cat(1,fullspontdffs2{:}); spontdists2 = cat(1,fullspontdists2{:});
%spontdffs = vertcat(spontdffs1,spontdffs2); spontdists = vertcat(spontdists1, spontdists2);
spontdffs = cat(1,fullspontdffs{:}); spontdists = cat(1,fullspontdists{:});
[binnums2, edges] = discretize(dists,bin_partition);
[binmeans2,binstds2, bincounts2] = grpstats(dffs2,binnums2,{'mean','std','numel'});
binstds2 = binstds2 ./ sqrt(bincounts2);

[spontbinnums, spontedges] = discretize(spontdists,bin_partition);
[spontbinnums2, spontedges2] = discretize(spontdists2,bin_partition);
[binmeansspont,binstdsspont,bincountsspont] = grpstats(spontdffs,spontbinnums,{'mean','std','numel'});
[binmeansspont2,binstdsspont2,bincountsspont2] = grpstats(spontdffs2,spontbinnums2,{'mean','std','numel'});
binstdsspont = binstdsspont ./ sqrt(bincountsspont);
ploterrbound =  binstdsspont;
empty_edges = setdiff(1:numel(edges),binnums2);
edges(empty_edges) = [];
hold on;
if plotspont
plot(spontedges(1:end-1), binmeansspont, 'x', 'MarkerEdgeColor', 'black', 'MarkerFaceColor','black' );
er = errorbar(edges, binmeansspont,ploterrbound, ploterrbound);
er.Color = 'black';
end
li = plot(edges,binmeans2, 'o', 'MarkerEdgeColor', 'g', 'MarkerFaceColor','g', 'MarkerSize',3);
li.Color = 'g';
er = errorbar(edges, binmeans2,binstds2, binstds2);
er.Color = 'g';
title('stim response vs stim distance');
axis(myaxes);

newxticks = 0:25: ceil(max(xticks)/(pix_um_ratio *25))*25;
xticks(newxticks * pix_um_ratio);
xticklabels(cellfun(@num2str,mat2cell(newxticks,1,ones(1,numel(newxticks))),'UniformOutput',false ));
%xticklabels(cellfun(@(x) num2str(floor(str2num(x) / pix_um_ratio)), xticklabels ,'UniformOutput',false));

grid on;
hold off;
axis square;

if adjscale
    yticks([0 .25 .5 .75 1]*adjmax);
    yticklabels([0 .25 .5 .75 1]);
end

% combine grpstats tables across experiments
grpstats_red = vertcat(grpstats_red{:});
grpstats_nonred = vertcat(grpstats_nonred{:});
grpstats_all = vertcat(grpstats_all{:});

dffsall = cat(1,alldffs{:}); distsall = cat(1,alldists{:});
[binnumsall, edgesall] = discretize(distsall,bin_partition);
[binmeansall,binstdsall, bincountsall] = grpstats(dffsall,binnumsall,{'mean','std','numel'});
binstdsall = binstdsall ./ sqrt(bincountsall);




% calculate p-values for chrimson +/- in bins

dffs = vertcat(dffs1,dffs2);
binnums = vertcat(binnums1, binnums2);
pvec_nonred_red=zeros(1,max(binnums));
pvec_nonred_spont = zeros(1,max(binnums));
pvec_red_spont =zeros(1,max(binnums));
pvec_nonred_allspont = zeros(1,max(binnums));
pvec_red_allspont = zeros(1,max(binnums));




for binnum=1:max(binnums)
   [~,pvec_nonred_red(binnum)] = ttest2(dffs1(binnums1==binnum),dffs2(binnums2==binnum));
   [~,pvec_nonred_spont(binnum)] = ttest2(spontdffs2(spontbinnums2==binnum),dffs2(binnums2==binnum));
   [~,pvec_red_spont(binnum)] = ttest2(dffs1(binnums1==binnum),spontdffs1(spontbinnums1==binnum));
   [~,pvec_red_allspont(binnum)] = ttest2(dffs1(binnums1==binnum),spontdffs(spontbinnums==binnum));
   [~,pvec_nonred_allspont(binnum)] = ttest2(dffs2(binnums2==binnum),spontdffs(spontbinnums==binnum));
end

sigtest_nonred_red = holm_bonferroni(pvec_nonred_red, .05);
sigtest_nonred_spont = holm_bonferroni(pvec_nonred_spont, .05);
sigtest_red_spont = holm_bonferroni(pvec_red_spont, .05);
sigtest_red_allspont = holm_bonferroni(pvec_red_allspont, .05);
sigtest_nonred_allspont = holm_bonferroni(pvec_nonred_allspont, .05);


binstats_red.count = bincounts1;
binstats_nonred.count = bincounts2;
binstats_all_spont.count = bincountsspont;
binstats_all.count = bincountsall;
binstats_red.mean = binmeans1;
binstats_nonred.mean = binmeans2;
binstats_all_spont.mean = binmeansspont;
binstats_all.mean = binmeansall;
binstats_nonred_spont.mean = binmeansspont2;
binstats_red_spont.mean = binmeansspont1;
binstats_red.sem = binstds1;
binstats_nonred.sem = binstds2;
binstats_nonred_spont.sem = binstdsspont2 ./ sqrt(bincountsspont2);
binstats_red_spont.sem = binstdsspont1 ./ sqrt(bincountsspont1);
binstats_all_spont.sem = binstdsspont;
binstats_all.sem = binstdsall;
bin_edges=  spontedges / pix_um_ratio;



binstats_red.cell_count = splitapply(@sum,grpstats_red.GroupCount > 0, grpstats_red.binnum);
binstats_nonred.cell_count = splitapply(@sum,grpstats_nonred.GroupCount > 0, grpstats_nonred.binnum);
binstats_all.cell_count = splitapply(@sum,grpstats_all.GroupCount > 0, grpstats_all.binnum);

[g,stim,binnum,expid] = findgroups(grpstats_red.stimcell,grpstats_red.binnum, grpstats_red.expid);
binstats_red.stim_count = splitapply(@sum,splitapply(@max,grpstats_red.GroupCount, g), binnum);
[g,stim,binnum,expid] = findgroups(grpstats_nonred.stimcell,grpstats_nonred.binnum, grpstats_nonred.expid);
binstats_nonred.stim_count = splitapply(@sum,splitapply(@max,grpstats_nonred.GroupCount, g), binnum);
[g,stim,binnum,expid] = findgroups(grpstats_all.stimcell,grpstats_all.binnum, grpstats_all.expid);
binstats_all.stim_count = splitapply(@sum,splitapply(@max,grpstats_all.GroupCount, g), binnum);


clearvars adj* spont* full* bin* -except binstats*;


%clear dffs1 dffs2 spontdffs spontdists dists fulldffs fulldists fullspontdffs fullspontdists

%% PCA on trace responses

alltraces = [];
k = Fcent_list.keys;
allstiminds = [prestimindsfull poststimindsfull];
for i=1:numel(k)
    kk = k{i};
    tmpmov = Fcent_list(kk);
    tmpmov = tmpmov(allstiminds,:,:);
    tmpmov = tmpmov ./ mean(tmpmov,1);
    tmpmov = reshape(tmpmov,[size(tmpmov,1) size(tmpmov,2)*size(tmpmov,3)]);
    alltraces = cat(2,alltraces,tmpmov);
end
alltraces = alltraces'; % rows are observations
R = corrcoef(alltraces);

[coeff,score,latent,tsquared,explained,mu] = pca(alltraces);

ncoeffs = 8;
%tic
%Y = tsne(score(1:5000,1:ncoeffs), 'NumDimensions',3, 'Distance','cosine');
%toc

%% combining local response maps between mice

blursigma = 5; % blurring filter
mapsize = 2*[height width];
mapcenter = round(mapsize/2);
finalzoom = 1.5;

celltype = 1; % set to [] for all cell types (except background)


clims_prob = [-.3,.3];
clims = [-.05 .05];
circle_cells = false; %if true, use approximating circles for cell shapes


spatial_map_excite = zeros(mapsize);
spatial_map_inhibit = zeros(mapsize);
spatial_map_dff = zeros(mapsize);
spatial_map_base = zeros(mapsize); %stores "denominator image" for map (local number of cells)



for iii=1:nmultifiles
    load(multi_files{iii},'stattab','isresponder','badcells', 'responders', 'xavg', 'yavg', 'maskinds');
    if numel(celltype)>0
        celltype_ind = stattab.celltype==celltype;
    else
        celltype_ind = stattab.celltype >= 0;
    end
    spatial_tab = grpstats(stattab(ismember(stattab.stimcell,find(isresponder)) & celltype_ind & ~ismember(stattab.cell,badcells) ,:),{'cell','stimcell','resptype'},{'mean'},'DataVars',{'dff','stimxdist','stimydist'});
    spatial_tab_temp.Properties.RowNames = {};

    
    for i=1:size(spatial_tab,1)
        icell = spatial_tab{i,'cell'};
        if circle_cells
            update_spot = spatial_spots(i);
        else
            [y,x] = ind2sub([height width], maskinds{icell});
            y = round(y - yavg(icell) - spatial_tab.mean_stimydist(i) + mapcenter(1));
            x = round(x - xavg(icell) - spatial_tab.mean_stimxdist(i) + mapcenter(2));
            update_spot = sub2ind(mapsize,y,x);
        end
        if spatial_tab.resptype(i)==1
            spatial_map_excite(update_spot) =  spatial_map_excite(update_spot) + spatial_tab.GroupCount(i);
        elseif spatial_tab.resptype(i)==-1
           spatial_map_inhibit(update_spot) =  spatial_map_inhibit(update_spot) + spatial_tab.GroupCount(i);
        end
        spatial_map_dff(update_spot) =  spatial_map_dff(update_spot) + spatial_tab.mean_dff(i)*spatial_tab.GroupCount(i);
        spatial_map_base(update_spot) = spatial_map_base(update_spot) + spatial_tab.GroupCount(i);
    end
end
    
   

    
spatial_map_excite =  imgaussfilt(spatial_map_excite, blursigma,'FilterSize',4*ceil(2*blursigma)+1.);
spatial_map_inhibit = imgaussfilt(spatial_map_inhibit, blursigma,'FilterSize',4*ceil(2*blursigma)+1.);
spatial_map_dff = imgaussfilt(spatial_map_dff, blursigma,'FilterSize',4*ceil(2*blursigma)+1.);
spatial_map_base = imgaussfilt(spatial_map_base, blursigma,'FilterSize',4*ceil(2*blursigma)+1.);

spatial_map_base(spatial_map_base == 0) = inf;
spatial_map_excite = spatial_map_excite ./ spatial_map_base;
spatial_map_inhibit = spatial_map_inhibit ./ spatial_map_base;
spatial_map_dff = spatial_map_dff ./ spatial_map_base ;


spatial_map_difference = spatial_map_excite - spatial_map_inhibit;


% center circle for plotting
th = 0:.1:2*pi+.1;
r = 20;
xunit = r * cos(th) + mapcenter(2);
yunit = r * sin(th) + mapcenter(1);

figure;
hold on;
imagesc(spatial_map_dff, clims); h=colorbar;colormap(mymap);
ylabel(h,'local mean DF/F0');
plot(xunit, yunit); title('spatial dff map');
xlim([mapcenter(2)-round(width/(2*finalzoom)),mapcenter(2)+round(width/(2*finalzoom))]); ylim([mapcenter(1)-round(height/(2*finalzoom)),mapcenter(1)+round(height/(2*finalzoom))]);

axis square
set(gca,'Layer','top');

newxticks = (-150:25:150); % in microns

xticks(mapcenter(2)+(newxticks * pix_um_ratio))
yticks(mapcenter(1)+(newxticks * pix_um_ratio))


xticklabels(cellfun(@num2str,mat2cell(newxticks,1,ones(1,numel(newxticks))),'UniformOutput',false ));
yticklabels(cellfun(@num2str,mat2cell(newxticks,1,ones(1,numel(newxticks))),'UniformOutput',false ));


c = redblue_log;
colormap(c);

xlabel('Displacement from stim site (microns)');
hold off;
%% PCA on cell-population trial DF/F vectors (scalar per cell)

%trialrecgrps = findgroups(stattab.trial, stattab.recid);
trialRespMatrix = unstack(stattab(ismember(stattab.stimcell,responders(4:8)),{'dff','recid','cell','trial','stimcell'}),'dff','cell');

[coeff,score,latent,tsquared,explained,mu] = pca(trialRespMatrix{:,4:end});

ncoeffs = 40;
tic
Y = tsne(score(:,1:ncoeffs), 'NumDimensions',3, 'Distance','cosine');
toc