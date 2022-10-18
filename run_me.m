%% step 0 - set up path 
clear all; close all; clc;
addpath(genpath(pwd));
gitDir = '/Users/pw1246/Documents/Github';
fsDir = '/Applications/freesurfer/7.2.0';
set_up(gitDir,fsDir)

subjid = 'sub-0201';
ses = 'ses-01';

% download example data here: 
% https://drive.google.com/file/d/1Ot5_QWl6whpB5qEXUEwQuzUJyis-oKJp/view?usp=sharing
% define path to the downloaded folder here:
bidsDir = '/Users/pw1246/Desktop/glm-Localizer-Tutorial-Data';

%% step 1 - load data and design matrix

load([bidsDir '/derivatives/data.mat']) % this is a 1x1 cell with size 302958 (vertices) by 300 (TR)

load dms  % this is a 1x1 cell with size 300 (TR) by 6 (conditions)

%% step 2 - take a look at the raw data and design matrix




%% step 3 - run kendrick's glm 

  myresults = GLMestimatemodel(dms,datafiles,1,1,'assume',[],0);
     
  
%% step 4 - take a look at R2, contrasts, and estimated betas

close all
figure(1);clf
bins = 0:1:100;
datatoplot = myresults.R2 ;
cmap0 = cmaplookup(bins,min(bins),max(bins),[],hot);

datatoplot(isnan(datatoplot)) = 0;

[rawimg,Lookup,rgbimg] = cvnlookup(subjid,1,datatoplot,[min(bins) max(bins)],cmap0,0,[],0,{'roiname',{'MT_exvivo'},'roicolor',{'w'},'drawrpoinames',0,'roiwidth',{5},'fontsize',20}); %MT_exvivo %Kastner2015

color = [0.5];
[r,c,t] = size(rgbimg);
[i j] = find(all(rgbimg == repmat(color,r,c,3),3));

for ii = 1: length(i)
    rgbimg(i(ii),j(ii),:) = ones(1,3);
end

a = imagesc(rgbimg); axis image tight;
axis off
hold on
% subplot(2,1,2)
plot(0,0);
colormap(cmap0);
hcb=colorbar('SouthOutside');
hcb.Ticks = [0:0.25:1];
% hcb.TickLabels = {'}
hcb.FontSize = 25
hcb.Label.String = 'R2%'
hcb.TickLength = 0.001;

title(subjid)

%%
figure(2); clf
betas = myresults.modelmd{2};

motion = nanmean(betas(:,[1]),2);
stationary = nanmean(betas(:,[2]),2);
 
C = [1 -1]'
contrast = C' * [motion stationary]';

alltcs = nanmean(cat(3,datafiles{:}),3);
predttcs = dms{1} * myresults.modelmd{2}';

datatoplot = contrast' .* double(myresults.R2>5);

datatoplot(datatoplot==0) = -50;
datatoplot(isnan(datatoplot)) = -50;

bins = -0.5:0.01:0.5
cmap0 = cmaplookup(bins,min(bins),max(bins),[],(cmapsign4));
[rawimg,Lookup,rgbimg] = cvnlookup(subjid,1,datatoplot,[min(bins) max(bins)],cmap0,min(bins),[],0,{'roiname',{'MT_exvivo'},'roicolor',{'w'},'drawrpoinames',0,'roiwidth',{5},'fontsize',20});

color = [0.5];
[r,c,t] = size(rgbimg);
[i j] = find(all(rgbimg == repmat(color,r,c,3),3));

for ii = 1: length(i)
    rgbimg(i(ii),j(ii),:) = ones(1,3);
end

a = imagesc(rgbimg); axis image tight;

    
set(gcf,'Position',[ 277         119        1141         898])
axis off
hold on
% subplot(2,1,2)
plot(0,0);
colormap(cmap0);
hcb=colorbar('SouthOutside');
hcb.Ticks = [0 1];
hcb.TickLabels = {'-0.5';'0.5'}
hcb.FontSize = 25
hcb.Label.String = 'Moving dots vs. stationary'
hcb.TickLength = 0.001;


%%
figure(3);clf

mymask = double(myresults.R2>min(maxk(myresults.R2,101)));
sum(mymask)

tcs = nanmean(cat(3,datafiles{:}),3);
ObsResp = nanmean(tcs(logical(mymask),:),1);
dc = nanmean(ObsResp)
ObsResp = 100 * (ObsResp - dc) / dc;

plot(ObsResp)
% tcs
hold on

stem(dms{1}(:,1))

legend({'Average timecourse from 100 most responsible voxels in MT';'Centrally moving dots predictior'})
legend box off
ylabel('%BOLD')
xlabel('TRs')
set(gca,'Fontsize',15)

 
%% step 5 - save the estimated responses into mgz files for future use
setenv('SUBJECTS_DIR',[bidsDir '/derivatives/freesurfer'])
conditions = {'central_moving';'central_stationary';'left_moving';'left_stationary';'right_moving';'right_stationary'}

resultsdir = sprintf('%s/derivatives/GLMdenoise/%s/%s/',bidsDir,subjid,ses);
fspth = fullfile(bidsDir, 'derivatives', 'freesurfer', subjid);
lcurv = read_curv(fullfile(bidsDir, 'derivatives', 'freesurfer', subjid, 'surf', 'lh.curv'));
rcurv = read_curv(fullfile(bidsDir, 'derivatives', 'freesurfer', subjid, 'surf', 'rh.curv'));
mgz = MRIread(fullfile(bidsDir, 'derivatives', 'freesurfer', subjid, 'mri', 'orig.mgz'));
leftidx  = 1:numel(lcurv);
rightidx = (1:numel(rcurv))+numel(lcurv);
%
for b = 1 : size(myresults.modelmd{2},2)
    
mgz.vol = myresults.modelmd{2}(leftidx,b);
MRIwrite(mgz, fullfile(resultsdir, sprintf('lh.%s.mgz',conditions{b})));
mgz.vol = myresults.modelmd{2}(rightidx,b);
MRIwrite(mgz, fullfile(resultsdir, sprintf('rh.%s.mgz',conditions{b})));
end

% findnan = isnan(myresults.R2(:));
% myresults.R2(findnan)= 0;
mgz.vol = double(myresults.R2(leftidx,:));
MRIwrite(mgz, fullfile(resultsdir, sprintf('lh.%s.mgz','vexpl_glm')));
mgz.vol = double(myresults.R2(rightidx,:));
MRIwrite(mgz, fullfile(resultsdir, sprintf('rh.%s.mgz','vexpl_glm')));


mgz.vol = double(myresults.R2(leftidx,:))>10;
MRIwrite(mgz, fullfile(resultsdir, sprintf('lh.%s.mgz','vexpl_mask')));
mgz.vol = double(myresults.R2(rightidx,:))>10;
MRIwrite(mgz, fullfile(resultsdir, sprintf('rh.%s.mgz','vexpl_mask')));


pairs =[[1 2];[3 4];[5 6]]


C = [1 -1]'

betas = myresults.modelmd{2};
for p = 1 : size(pairs,1)
    
motion = nanmean(betas(:,[pairs(p,1)]),2);
stationary = nanmean(betas(:,[pairs(p,2)]),2);    
contrast = (C' * [motion stationary]')';

mgz.vol = double(contrast(leftidx,:));
MRIwrite(mgz, fullfile(resultsdir, sprintf('lh.%s_vs_%s.mgz',conditions{pairs(p,1)},conditions{pairs(p,2)})))
mgz.vol = double(contrast(rightidx,:));
MRIwrite(mgz, fullfile(resultsdir, sprintf('rh.%s_vs_%s.mgz',conditions{pairs(p,1)},conditions{pairs(p,2)})))
end

%% step 6 - draw some ROIs


roilabels = {'hMT';'MST';}
rng    = [0 2]
resultsdir = sprintf('%s/derivatives/GLMdenoise/%s/%s/',bidsDir,subjid,ses);
mgznames = {'central_moving_vs_central_stationary';'left_moving_vs_left_stationary';'right_moving_vs_right_stationary'}; 
for zz = 1:3
cmap = cmapsign4;
bins = [-0.5:0.1:0.5];
cmaps{zz} = cmaplookup(bins,min(bins),max(bins),[],cmap);
crngs{zz} = [min(bins) max(bins)];
threshs{zz} = min(bins);
end
drawMyRois



