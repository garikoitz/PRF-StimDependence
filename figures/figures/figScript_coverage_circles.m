
%% Plot the FOV a la Kay et al. figures:
% Randomly sample 100 voxels from an ROI
% Plot the centers and draw the circle for the half max

clear all; close all; clc; 
bookKeeping; 

%% modify here
% subjects to include
list_subInds = 	[13:20];

% threshold and other visualization parameters
vfc = ff_vfcDefault; 
vfc.cothresh = 0.45

% rois to do this for
list_roiNames = {
    'WangAtlas_V1v_left'
    'WangAtlas_hV4_left'
    'WangAtlas_IPS0_left'
    'lVOTRC'
    };

% modify this if roiIndividual
dotColor = [1 0 0];

% ret model dt and name
list_dtNames = {
    'WordLarge'
    'WordSmall'
    };
list_rmNames = {
%     'retModel-Words-css.mat'
%     'retModel-Checkers-css.mat'
    'retModel-WordLarge-css.mat'
    'retModel-WordSmall-css.mat'
    }; 

% number of voxels to sample
numVoxels = 100; 

% whether to color each pRF center red, or to differentiate
colorRed = false; 
if ~colorRed
    voxelColors = jetCmap(0,numVoxels);
end

%% initialize
numRois = length(list_roiNames);
numRms = length(list_rmNames);
numSubs = length(list_subInds);

%% get the data once
rmroiCellUnthresh = ff_rmroiCell(list_subInds, list_roiNames, list_dtNames, list_rmNames);

%% Threshold and get identical voxels for each subject
% In comparing ret models, the collection of voxels may not be the same
% because of the thresholding. In this cell we redefine the rmroi
rmroiCell = cell(size(rmroiCellUnthresh));

for jj = 1:numRois
    for ii = 1:numSubs        
        % get identical voxels for each subject's roi over all ret models
        D = rmroiCellUnthresh(ii,jj,:);
        rmroiCell(ii,jj,:) = ff_rmroiGetSameVoxels(D, vfc);        
    end
end


close all;

for jj = 1:numRois
    
    roiName = list_roiNames{jj};
    
    for kk = 1:numRms
       figure;
       dtName = list_dtNames{kk};
        
       %% linearize the sigma and ecc data
       rmroiSubjects = rmroiCell(:,jj,kk);       
       ldata_ecc = ff_rmroiLinearize(rmroiSubjects, 'ecc');
       ldata_sigmaEff = ff_rmroiLinearize(rmroiSubjects, 'sigma');
       ldata_sigmaMaj = ff_rmroiLinearize(rmroiSubjects, 'sigma1');
       ldata_exponent = ff_rmroiLinearize(rmroiSubjects, 'exponent');
       
       ldata_x0 = ff_rmroiLinearize(rmroiSubjects, 'x0');
       ldata_y0 = ff_rmroiLinearize(rmroiSubjects,'y0');
       
       %% randomly sample
       % 100 indices
       vecInd = [1:length(ldata_ecc)]; 
       randInd = datasample(vecInd,numVoxels, 'replace', false); 
       
       x100 = ldata_x0(randInd); 
       y100 = ldata_y0(randInd);
       sigma100 = ldata_sigmaEff(randInd);
       
       %% plotting 
       axis([-vfc.fieldRange vfc.fieldRange -vfc.fieldRange vfc.fieldRange])
       ff_polarPlot(vfc);
       
       % the circles
       % formula
       hold on; 
       viscircles([x100' y100'], sigma100, 'Color', [.1 .1 .1], 'Linewidth',1)
       
       % the centers
       if colorRed
           scatter(x100,y100, 'r','filled')
       else
           scatter(x100, y100, [], voxelColors, 'filled');
       end
       
       
       % titles
       titleName = {
        '100 randomly selected pRFs'
        [roiName '. ' dtName]
        };
       title(titleName, 'fontweight', 'bold')
       
    end
    
end