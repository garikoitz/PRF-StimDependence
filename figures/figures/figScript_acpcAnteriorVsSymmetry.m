%% Plot symmetry: (ipsi-contra) as a function voxel atneriorness

% clear all; close all; clc; 
% bookKeeping; mrvCleanWorkspace

%% modify here

list_subInds = 1:20;
list_paths = cr.bk.list_sessionRet; 

roiName = 'lVOTRC-threshByWordModel';

dtName = 'Words';
rmName = 'retModel-Words-css.mat';

% which hemifield. 'left' | 'right'
wHemifield = 'left'; 

vfc = ff_vfcDefault; 
vfc.cothresh = 0; 
contourLevel = 0.5; 

%% do the things

figure; hold on; 

for ii = list_subInds
    
    dirAnatomy = cr.bk.list_anatomy{ii};
    dirVista = list_paths{ii};
    chdir(dirVista);
    subColor = cr.bk.list_colorsPerSub(ii,:);
    subInitials = cr.bk.list_sub{ii};
    
    %% init gray, load roi
    vw = initHiddenGray; 
    roiPath = fullfile(dirAnatomy, 'ROIs', roiName);
    vw = loadROI(vw, roiPath, [],[],1,0); 
    
    % the coordinates of the ROI
    % different from voxels ... so many more
    % for now, just run with this
    roiCoords = vw.ROIs.coords; 
    
    % number of coordinates in this ROI
    numCoords = size(roiCoords,2);
    
    %% load the ret model and get the rmroi
    rmPath = fullfile(dirVista, 'Gray', dtName, rmName); 
    vw = rmSelect(vw,1,rmPath);
    vw = rmLoadDefault(vw);
    rmroi = rmGetParamsFromROI(vw);
        
    %% loop over roi coords
    for nn = 1:numCoords
        
        % make a point ROI in the view ... this should automatically select
        coord = roiCoords(:,nn); 
        vw = makePointROI(vw, coord, 1); 
        
        %% get acpc coordinates
        % if skipTalFlag is 0, we use talairach coordinages
        % if skipTalFlag is 1, we use the spatial norm
        skipTalFlag = 1; % default is 0
        snc = vol2talairachVolume(vw,skipTalFlag);
        
        % the y coordinate connects "ac" to "pc"
        % store the y coordinates
        % the more negative the number, the more posterior it is
        anteriorness = snc(2); 
        
        % delete the roi so we don't clog up the view
        vw = deleteROI(vw,2);
        
        
        %% plotting
        x0 = rmroi.x0(nn);
        y0 = rmroi.y0(nn);
        sig = rmroi.sigma(nn);
        
        leftArea = ff_hemifieldArea(x0,y0,sig, contourLevel,'left', vfc); 
        rightArea = ff_hemifieldArea(x0,y0, sig, contourLevel, 'right', vfc); 
        areaDifference = leftArea - rightArea; 
        
        if ~isnan(sig)
            
            plot(anteriorness,areaDifference,'o', ...
                'MarkerFaceColor', subColor, 'MarkerEdgeColor', [0 0 0])
        end
      
    end
 
end

xlabel('MNI coordinate. posterior --> anterior')
ylabel(['Left area - Right area (deg2)'])
titleName = {
    'Anteriorness and symmetry'
    [roiName]
    };
title(titleName, 'fontweight','bold')
grid on; 

% save
ff_dropboxSave; 
