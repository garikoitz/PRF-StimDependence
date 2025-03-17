%% INIT
tbUse PRF-StimDependence;

close all; clear all;
mrvCleanWorkspace;
cr         = struct();
cr.codeDir = sdRP;

% WHERE THE NEW DATA IS
cr.dirs.BASE     = '/acorn/data/neuro/gari/PRF-StimDependence';
cr.dirs.DATA     = fullfile(cr.dirs.BASE,'DATA');
cr.dirs.ANALYSIS = fullfile(cr.dirs.BASE,'ANALYSIS');
cr.dirs.ORG      = fullfile(cr.codeDir,'DATA','ANALYSIS','matlabfiles','organization');
cr.dirs.DEF      = fullfile(cr.codeDir,'DATA','ANALYSIS','matlabfiles','defineProjectDefaults');
cr.dirs.FIG     = fullfile('~/toolboxes/PRF-StimDependence','DATA', 'figures');

cr.dirs.FIGPNG  = fullfile(cr.dirs.FIG,'png');
cr.dirs.FIGSVG  = fullfile(cr.dirs.FIG,'svg');

cr.bk = bookKeeping(cr);


%% GENERIC
list_roiNames = {'WangAtlas_V1d_left'
                 'WangAtlas_V2d_left'
                 'WangAtlas_V3d_left'
                 'WangAtlas_V1v_left'
                 'WangAtlas_V2v_left'
                 'WangAtlas_V3v_left'
                 'WangAtlas_hV4_left'
                 'WangAtlas_VO1_left'
                 'WangAtlas_V3A_left'
                 'WangAtlas_IPS0_left'
                 'WangAtlas_IPS1_left'
                 'WangAtlas_V1d_right'
                 'WangAtlas_V2d_right'
                 'WangAtlas_V3d_right'
                 'WangAtlas_V1v_right'
                 'WangAtlas_V2v_right'
                 'WangAtlas_V3v_right'
                 'WangAtlas_hV4_right'
                 'WangAtlas_VO1_right'
                 'WangAtlas_V3A_right'
                 'WangAtlas_IPS0_right'
                 'WangAtlas_IPS1_right'};

new_roi_ind = [1,2,3,4,5,6,7,8,9,10,11]; % left only
% new_roi_ind = [12:length(list_roiNames)]; % right only
% new_roi_ind = [1:length(list_roiNames)]; % all 
new_list_roiNames = list_roiNames(new_roi_ind); 
               
% Read the generic params for coverage for all subjects
cr.defaults.covfig.vfc = ff_vfcDefault();
cr.defaults.covfig.vfc.list_roiNames = new_list_roiNames;

% CNI
CNI_list_subInds = [1:20];
CNI_data_type_ind = [2, 1]; % to have words in x and CB in y
CNI_rmroiFname='rmroicell_subInds-1to20_dtNames-cb-w-ff_fits-new_LeftRightROIs_2023.mat';
CNI_fpath = fullfile(sdRP,'DATA',CNI_rmroiFname);
CNI_rmroiCell = load(CNI_fpath,'rmroiCell').rmroiCell;
CNI_list_dtNames     = {'Checkers','Words','FalseFont'};
CNI_list_rmDescripts = CNI_list_dtNames;
% make it WORD VS CB
CNI_data_types = CNI_list_dtNames(CNI_data_type_ind);
new_CNI_rmroiCell = CNI_rmroiCell(:, new_roi_ind, CNI_data_type_ind);
cr.defaults.covfig.vfc.list_dtNames = CNI_list_dtNames;
varExplained=0.2;
fieldrange = 15;
[CNI_R,CNI_C_data,CNI_cr]=crThreshGetSameVoxel( cr,...
                                    new_CNI_rmroiCell,...
                                    CNI_list_subInds,...
                                    new_list_roiNames,...
                                    'cothres', varExplained,...
                                    'fieldrange', fieldrange, ...
                                    'show_summary', false);


% Create distribution comparisons plots
save_fig = true;
show_fig = 'off';
filter_value = 5; % deg of visual field
fname = ['distr_comp_CNI_RIGHT_' CNI_data_types{2} '-'  CNI_data_types{1}];
fpath = '~/toolboxes/PRF-StimDependence/DATA/figures/png/distr';
path_fname = fullfile(fpath, [fname '.png']);
sd_create_distr_plot(CNI_R, ...
                     CNI_data_types, ...
                     new_list_roiNames, ...
                     varExplained, ...
                     'save_fig' , save_fig, ...
                     'path_fname' , path_fname, ...
                     'filter_value' , filter_value, ...
                     'show_fig' , show_fig);

% HEB
HEB_list_subInds  = [31:36 38:44];
HEB_data_type_ind = [1, 3]; % to have words in x and CB in y
HEB_rmroiFname='rmroicell_subInds-31to36-38to44_dtNames-ALL-LeftRight_fits-new_2023.mat';
HEB_fpath = fullfile(sdRP,'DATA',HEB_rmroiFname);
HEB_rmroiCell = load(HEB_fpath,'rmroiCell').rmroiCell;
HEB_list_dtNames     = {'Words_English','Words_Hebrew','Checkers'};
HEB_list_rmDescripts = HEB_list_dtNames;
% WORD VS CB
HEB_data_types = HEB_list_dtNames(HEB_data_type_ind);
new_HEB_rmroiCell = HEB_rmroiCell(:,new_roi_ind,CNI_data_type_ind);

cr.defaults.covfig.vfc.list_dtNames = HEB_list_dtNames;
varExplained=0.05;
fieldrange = 7;
[HEB_R,HEB_C_data,HEB_cr]=crThreshGetSameVoxel( cr,...
                                    new_HEB_rmroiCell,...
                                    HEB_list_subInds,...
                                    new_list_roiNames,...
                                    'cothres', varExplained,...
                                    'fieldrange', fieldrange, ...
                                    'show_summary', false);

% Create distribution comparisons plots
save_fig = true;
show_fig = 'off';
filter_value = 2.34; % deg of visual field
fname = ['distr_comp_HEB_LEFT_' HEB_data_types{2} '-filterdeg-' num2str(filter_value) '_'  HEB_data_types{1}];
fpath = '~/toolboxes/PRF-StimDependence/DATA/figures/png/distr';
path_fname = fullfile(fpath, [fname '.png']);
sd_create_distr_plot(HEB_R, ...
                     HEB_data_types, ...
                     new_list_roiNames, ...
                     varExplained, ...
                     'save_fig' , save_fig, ...
                     'path_fname' , path_fname, ...
                     'filter_value' , filter_value, ...
                     'show_fig' , show_fig);

% CHECK AND NOTE HERE THE DIFFERENCES FOR ATTENTION (Kendrick, Kalanit, check brian jon paper review, ...)


%{                                   

%% See time series of a subject that has higher CB and lower 
% Sub:1 is heb_pilot_09
%{ 
% Sub:1, roi:WangAtlas_hV4_left, ecc_diff: 0.94, coord:125 175 68, index:192395 
subj = '09';
roi = 'hV4';
ind = 192395;
%}
%{ 
% Sub:5, roi:WangAtlas_hV4_left, ecc_diff: 9.1, coord:122 171 74, index:147830
subj = '13'; % 4 + 9
roi = 'hV4';
ind = 147830;
%}
%{ 
% Sub:8, roi:WangAtlas_hV4_left, ecc_diff: 2.3, coord:135 190 67, index:230801 
subj = '16'; % 7 + 9
roi = 'hV4';
ind = 230801;
%}


plot_params = struct();
plot_params.p2_ret_data = ['/acorn/data/neuro/gari/PRF-StimDependence/DATA/heb_pilot' subj '/RetAndHebrewLoc/Gray'];
plot_params.generic = 'TSeries/Scan1/tSeries1.mat';
plot_params.ind = ind;
plot_params.what_data_types =  {'Words_English' ,'Words_Hebrew','Checkers'};
plot_params.rmNames = {'retModel-Words_English-css-fFit-fFit.mat',...
           'retModel-Words_Hebrew-css-fFit-fFit.mat',...
           'retModel-Checkers-css-fFit-fFit.mat'};


fname = ['heb_sub' subj '_' roi '_left_ind-' num2str(ind)];
fpath = '~/toolboxes/PRF-StimDependence/DATA/figures/png';
path_fname = fullfile(fpath, [fname '.png']);
plot_time_series(plot_params, path_fname)

%}

%{
 print legend with the data type and R2 rmse
 bootstrap the data to estiamte how confident we are on the model fit. 
  - we are minig


  ask david to send the decimation experiment

%}