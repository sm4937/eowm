function [RFs,anti_RFs] = draw_RFs(PRFparams,ROI,type)
%draw_RFs Draw voxel-specific receptive fields based on the PRF params
% provided
%  Take means from PRF params, stds, draw Gaussian receptive fields, as
%  they are assumed to be in mrVista fitting algorithm

nvoxels = sum(ROI.data(:));

% Note: To find PRF params for each subject, look for file names
% resembling:
% PRFparams = niftiread('/System/Volumes/Data/d/DATA/data/vRF_tcs/CC/RF1/CC_RF1_vista/RF_ss5_25mm-fFit.nii.gz');

% Params is 8D, and each dimension is:
% pol : polar angle
% ve : variance explained
% ecc : eccentricity
% sigmamajor : sigma of the gaussian describing the RF
% exponent : the exponent of the gaussian RF
% x0 : center of the PRF in x coords
% y0 : center of the PRF in y coords, flipped from what you would think
% i.e. negative number are above the horizontal meridian
% b : amplitude

voxel_RF_centers = rad2deg(PRFparams(:,:,:,1)); 
voxel_RF_centers = voxel_RF_centers(ROI.data>0);
voxel_RF_sigmas = PRFparams(:,4);

outRF_means = NaN(nvoxels,1);

antiRF = (1:360)-180;
antiRF(antiRF<1) = antiRF(antiRF<1) + 360;
for ii = 1:length(voxel_RF_centers)
    idx = ceil(voxel_RF_centers(ii));
    idx(idx > 360) = idx(idx>360)-360; idx(idx<1) = idx(idx<1)+360;
    outRF_means(ii,:) = antiRF(idx);
end

if strcmp(type,'box')
    % Do a box car RF function with a set number of degrees on either side

    RFs = zeros(nvoxels,360); anti_RFs = zeros(nvoxels,360);
    surround = 15;

    for voxel = 1:nvoxels
        temp = false(1,360);
        indices = round(voxel_RF_centers(voxel))-surround:round(voxel_RF_centers(voxel))+surround;
        indices(indices<1) = indices(indices<1)+360;
        indices(indices>360) = indices(indices>360) - 360;
        %make circular RFs, not just abrupt endings at edges
        temp(:,indices) = true;
        % assign RFs to specific voxel
        RFs(voxel,:) = temp;

        temp = false(1,360);
        indices = round(outRF_means(voxel))-surround:round(outRF_means(voxel))+surround;
        indices(indices<1) = indices(indices<1)+360;
        indices(indices>360) = indices(indices>360) - 360;
        temp(:,indices) = true;
        anti_RFs(voxel,:) = temp;
    end
    
else strcmp(type,'bump') 
    % Do a smooth Gaussian bump, weighting contribution of RF by proximity,
    % width of RF
    
    RFs = normpdf(
end
    
end
    
end

