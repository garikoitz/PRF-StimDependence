function fg = IndivFigDefaults
%INDIVFIGDEFAULTS Parameter default to run 
%   I think this won't be very useful, as it is project dependent. 
% When calling the function figFunction_coverage_individual() it will require
% some defaults, and this file will be passed, but the important thing is that
% every time this package is run, we will need to provide the proper defaults,
% and those will be project specific. 

% For now leave all the options here. This is generating the correct plots for
% reproducing them in the paper. 

fg.titleDescript = 'FOV';

% vfc threshold
fg.vfc = ff_vfcDefault_Hebrew; 
fg.vfc.cothresh = 0.2; 
fg.vfc.cmap = 'hot';
fg.vfc.addCenters = true;
fg.vfc.contourPlot = true;

% subjects
% list_subInds =  [31:36 38:44];
fg.list_subInds =  [31];

% session
% list_sessionHebrewRet, list_sessionRet
fg.list_path = list_sessionRet; 

% roi
% lh_VWFA_rl
% lh_VWFA_fullField_rl
fg.list_roiNames = {
%     'WangAtlas_V1'
%     'WangAtlas_V2'
%     'WangAtlas_V3'
%     'WangAtlas_hV4'
%     'LV2v_rl'
    'lVOTRC'
%     'lVOTRC-threshBy-Words_EnglishOrWords_English-co0p05'
%     'lVOTRC'
%     'rVOTRC'
%     'left_VWFA'
%     'lh_ventral_3_rl'
%     'combined_VWFA_rl'
%     'right_VWFA_rl'
%     'lh_VWFA_rl'
%     'lh_VWFA_fullField_rl'
%     'LV1_rl'
%     'LV2v_rl'
%     'LV3v_rl'
%     'rVOTRC'
    };

% dt and rm names
fg.list_dtNames = {
    'Words_Hebrew'
    'Words_English'
%     'Words_English'
%     'Words_English'
%     'Checkers'
%     'Words'
%     'Words1'
%     'Words2'
%     'Words_scale1mu0sig1'
%     'Words_scale1mu0sig0p5'
%     'Words_scale0p5mu0sig0'
%     'Checkers'
    };
fg.list_rmNames = {
    'retModel-Words_Hebrew-css.mat'
    'retModel-Words_English-css.mat'
%     'retModel-Words_English-css-4p5rad.mat'
%     'retModel-Words_English-css-4p5rad.mat'
%     'retModel-Checkers-css.mat'
%     'retModel-Words-css.mat'
%     'retModel-Words1-css.mat'
%     'retModel-Words2-css.mat'
%     'retModel-Words_scale1mu0sig1-css-left_VWFA_rl.mat'
%     'retModel-Words_scale1mu0sig0p5-css-left_VWFA_rl.mat'
%     'retModel-Words_scale0p5mu0sig0-css-left_VWFA_rl.mat'
%     'retModel-Checkers-css.mat'
    };

fg.list_rmDescripts = {
    'Words_Hebrew'
    'Words_English'
    };
end

