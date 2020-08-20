function vfc = ff_vfcDefault_Hebrew()

% vfc threshold
vfc.prf_size        = true; 
vfc.fieldRange      = 7;
vfc.method          = 'max';         
vfc.newfig          = true;                      
vfc.nboot           = 50;                          
vfc.normalizeRange  = true;              
vfc.smoothSigma     = false;                
vfc.cothresh        = 0.2;         
vfc.eccthresh       = [0 7];
vfc.sigmaEffthresh  = [0 7];   % sigma effect (sigmaMajor/sqrt(exponent))
vfc.sigmaMajthresh  = [0 14];  % sigma major 
vfc.nSamples        = 128;            
vfc.meanThresh      = 0;
vfc.weight          = 'fixed';  
vfc.weightBeta      = false;
vfc.cmap            = 'jet';						
vfc.clipn           = 'fixed';                    
vfc.threshByCoh     = true;                
vfc.addCenters      = false;                 
vfc.verbose         = prefsVerboseCheck;
vfc.dualVEthresh    = 0;
vfc.backgroundColor = [.9 .9 .9]; 
vfc.ellipsePlot     = false; 
vfc.ellipseLevel    = 0.5;
vfc.ellipseColor    = [1 0 0];
vfc.contourPlot     = true; 
vfc.contourLevel    = 0.5; 
vfc.contourColor    = [0 0 0];
vfc.tickLabel       = false; 
vfc.contourBootstrap= false; 
vfc.cothreshceil    = 1;        % don't get voxels higher than this
vfc.gridColor       = [.6 .6 .6];



% rois we want to look at
vfc.list_roiNames = {
%     'LV2v_rl'
%     'WangAtlas_V1d_left'
%     'WangAtlas_V2v_left'
%     'WangAtlas_V3v_left'
%     'WangAtlas_hV4_left'
%     'WangAtlas_VO1_left'
%     'lVOTRC'
%     'WangAtlas_V1v.mat'
%     'WangAtlas_V2v.mat'
%     'WangAtlas_V3v.mat'
%     'WangAtlas_hV4.mat'
%     'WangAtlas_IPS1'
%     'WangAtlas_IPS2'
%     'WangAtlas_IPS3'
%     'WangAtlas_IPS4'
%     'WangAtlas_IPS5'
%      'cVOTRC'
     'lVOTRC'
     'rVOTRC'
%     'Cohen2002VWFA_8mm'
%     'VOTRC'
%     'left_FFAFace_rl' 
%     'right_FFAFace_rl'
%     'LhV4_rl-threshByWordModel'
%     'LhV4_rl-threshByWordModel'
%     'RV1_rl'
%     'RV2v_rl'
%     'RV3v_rl'
%     'combined_VWFA_rl'
%     'right_VWFA_rl'
%     'ch_PPA_Place_rl'
%     'rh_PPA_Place_rl'
%     'lh_ventral_3_rl'
%     'lh_pFusFace_rl'
%     'lh_mFusFace_rl'
%     'lVOTRC'
%     'rVOTRC-threshByCheckerModel'
%     'rVOTRC-threshByWordModel'
%     'LV1_rl'
    };

% data types we want to look at
vfc.list_dtNames = {
    'Checkers'
    'Words_Hebrew'
    'Words_English'
%     'Words'
%     'Checkers'   
%     'WordLarge'
%     'WordSmall'
%     'FaceSmall'
%     'FaceLarge'
%     'Checkers'
%     'Words1'
%     'Words2'
%     'Words_scale1mu0sig1p5'
    };

% names of the rm in each dt
vfc.list_rmNames = {
    'retModel-Checkers-css.mat'
    'retModel-Words_Hebrew-css.mat'
    'retModel-Words_English-css.mat'
%     'retModel-Checkers-css.mat'
%     'retModel-WordLarge-css.mat'
%     'retModel-WordSmall-css.mat';
%     'retModel-FaceSmall-css.mat';
%     'retModel-FaceLarge-css.mat';
%     'retModel-Checkers-css.mat'
%     'retModel-Words1-css.mat'
%     'retModel-Words2-css.mat'
%      'retModel-Words_scale1mu0sig1p5-css-rVOTRC.mat'
    };






end