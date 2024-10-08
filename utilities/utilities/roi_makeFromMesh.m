%% script to define a new roi based on functional activity on mesh
% run this script after drawing and filling an roi
% rl, 08/2014


%% modify

% which view variable
% 1: vw
% 2: VOLUME{end}
% 3: hG
wView = 2;

newroi.color    = 'b';
newroi.name     = 'rVOTRC'; % 'lh_VWFA_fullField_WordVFaceScrambled_rl' 'GLM_WordVFace_Scrambled'
newroi.comment  = '';
restrictToFunc  = 1;  % 0 for visual field maps, 1 for categories
saveWhere       = 0; % 1 = lmocal, 0 = shared
% 'lh_VWFA_rl'          : on the inferior temporal sulcus (iTS), posterior of
%       the mid-fusiform suclus. also sometimes on posterior fusiform
%       gyrus. the iTS kind of curves upwards and is L-shaped
% 'lh_OWFA_rl'          : anything stemming from confluent fovea
% 'lh_WordsVentral_rl' : combine VWFA and OWFA    
% 'lh_WordVAll_rl'      : All voxels that are activated when making this

% 'lh_mFus_Face_rl'     : mid fusiform. medial or on the OTS.
% 'lh_pFus_Face_rl'     : posterior fusiform. anterior or on the pTCS.
% 'lh_iOG_Face_rl'      : posterior of the pTCS
% 'lh_FacesVentral_rl'  : everything on the ventral surface

%% no need to modify
% view with a loaded mesh
% 1: vw
% 2: VOLUME{end}
% 3: hG
switch  wView
    case 1
       theView = vw;  
    case 2
       theView = VOLUME{end}; 
    case 3
       theView = hG; 
end


% get the roi from the mesh
% vw = meshROI2Volume(vw, [mapMethod=3]), where method 3 means grow from 
% layer 1 to get an roi that spans all layers  
theView = meshROI2Volume(theView, 3); 

% whether or not to restrict roi to functional acitivity
if restrictToFunc
    theView = restrictROIfromMenu(theView); 
end

% grab selected roi

roi = viewGet(theView, 'curRoi'); 

%% perform ROI a not b
% roi is last one you picked
roiA = theView.ROIs(end).name;

% all other rois you don't want
roiB={theView.ROIs(1:end-1).name}; 

% make the roi 
theView = ff_ROIanotb(theView, roiA, roiB, newroi.name, newroi.color); 


%% save roi in local directory
saveROI(theView, 'selected', 0)

% refresh screen
theView = refreshScreen(theView); 
% refresh mesh
theView = meshColorOverlay(theView); 


% %% plot the coverage
% coverage_plot