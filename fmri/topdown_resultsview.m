%% Preview group-level Top-Down analysis results
% Code written by Shanshan Li May 2023

results = niftiRead('/System/Volumes/Data/d/DATB/datb/eowm_SM/TopDown_files/group_t&p_against0.nii');
stats = results.data;
tmapn = stats(:,:,:,1);
pmapn = stats(:,:,:,2);

thresholdtmap = tmapn;
thresholdtmap(pmapn > 0.05) = NaN;
figure;
subplot(1,2,1)
vol3d('cdata', tmapn, 'texture', '3D');
view(3);
axis tight;
daspect([1 1 1]);
colormap jet;
colorbar;
title('tmap');
subplot(1,2,2)
vol3d('cdata', thresholdtmap, 'texture', '3D');
view(3);
axis tight;
daspect([1 1 1]);
colormap jet;
colorbar;
title('threshold, p < 0.05');