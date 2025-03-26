clear all; 
% install tbUse
addpath(genpath('/acorn/data/neuro/soft/toolboxes/PRF-StimDependence'));
addpath(genpath('/acorn/data/neuro/soft/toolboxes/vistasoft'));
addpath(genpath('/acorn/data/neuro/soft/toolboxes/scitran'));


%% Time series test
cr         = struct();
cr.codeDir = sdRP;
cr.dirs.BASE     = '/acorn/data/neuro/gari/PRF-StimDependence';
cr.dirs.DATA     = fullfile(cr.dirs.BASE,'DATA');
cr.dirs.ANALYSIS = fullfile(cr.dirs.BASE,'ANALYSIS');
cr.dirs.ORG      = fullfile(cr.codeDir,'DATA','ANALYSIS','matlabfiles','organization');
cr.dirs.DEF      = fullfile(cr.codeDir,'DATA','ANALYSIS','matlabfiles','defineProjectDefaults');
DATA_dir = cr.dirs.DATA;
cr.bk = bookKeeping(cr);
cd(DATA_dir);

% normalize function from matlab is being shadowed by another one with same name in here
% rmpath(fullfile(fileparts(sdRP),'vistasoft/mrBOLD/Utilities/'));

subinds = [31:36 38:44];
% list_subInds        = [31:36 38:44];  % Hebrew
do_all_sub_plot = true; % if true, plot all subjects in one figure
do_indiv_sub_plot = false;
roi = 'V2v';
hemi = 'left';
% res_end = '-fFit-fFit'; % check if this was analuzed with the right fov
% res_end = '-fFit';
res_end = '';

%{
    HEBREW
    'heb_pilot09' = subind 31 = Subject 01
    'heb_pilot10' = subind 32 = Subject 02
    'heb_pilot11' = subind 33 = Subject 03
    'heb_pilot12' = subind 34 = Subject 04
    'heb_pilot13' = subind 35 = Subject 05
    'heb_pilot14' = subind 36 = Subject 06
%} 
% subname = 'heb_pilot13'; ind 35
% anat_subname = 'Dael';
% scan_name = 'RetAndHebrewLoc';
if do_all_sub_plot
    all_CB_ecc = []; all_WE_ecc = []; all_WH_ecc = [];
    all_CBWE_corr = []; all_CBWH_corr = []; all_WEWH_corr = []; 
end

% subject 5, VO1, index 156949, from functions/compare_time_series
% my_voxel = 156949; > calculate from this selection, not from compare_time_series
% save image dir
image_dir = fullfile(sdRP,'DATA','figures','png');

% FILTERS
% To filter too small VE
VE = .05;
% To eliminate de fovea
min_ecc = 0;
% To filter too small size
min_size = 0;


for subind = subinds
    mrvCleanWorkspace;
    scan_dir = cr.bk.list_sessionRet{subind};
    gray_dir = fullfile(scan_dir,'Gray');
    anat_dir = cr.bk.list_anatomy{subind};
    [~,b]=fileparts(fileparts(scan_dir));
    [~,c]=fileparts(anat_dir);
    subname = [num2str(subind) '-' b '-' c];
    sprintf('\n\nProcessing %s, ROI: %s_%s, result: %s\n', ...
            subname, hemi, roi, res_end)

    % Load ROIs
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
    roiPath = fullfile(anat_dir,'ROIs',['WangAtlas_' roi '_' hemi '.mat']);
    AA      = load(roiPath);
    % Load the vw
    cd(scan_dir)
    setpref('VISTA', 'fileFormat', 'nifti')
    load(fullfile(scan_dir,'mrSESSION.mat'))
    t1_file_path = fullfile(anat_dir,'t1.nii.gz');
    setVAnatomyPath(t1_file_path);
    vw = initHiddenGray;
    vw = rmLoadDefault(vw);
    vw    = loadROI(vw, roiPath, [],[],1,0);
    assert(isequal(AA.ROI.coords, vw.ROIs(1).coords));

    ROIindices   = viewGet(vw, 'roiIndices');

    % Read results
    CBresults = load(fullfile(gray_dir, ...
                    'Checkers',['retModel-Checkers-css' res_end '.mat']));
    [~,CBecc] = cart2pol(CBresults.model{1}.x0, CBresults.model{1}.y0); 
    CBr2 = 1 - CBresults.model{1}.rss ./ CBresults.model{1}.rawrss;
    CBsize = CBresults.model{1}.sigma.major;

    WEresults = load(fullfile(gray_dir,...
                    'Words_English','retModel-Words_English-css-fFit-fFit.mat'));
    [~,WEecc] = cart2pol(WEresults.model{1}.x0, WEresults.model{1}.y0); 
    WEr2 = 1 - WEresults.model{1}.rss ./ WEresults.model{1}.rawrss;
    WEsize = WEresults.model{1}.sigma.major;

    WHresults = load(fullfile(gray_dir,...
                    'Words_Hebrew','retModel-Words_Hebrew-css-fFit-fFit.mat'));
    [~,WHecc] = cart2pol(WHresults.model{1}.x0, WHresults.model{1}.y0); 
    WHr2 = 1 - WHresults.model{1}.rss ./ WHresults.model{1}.rawrss;
    WHsize = WHresults.model{1}.sigma.major;


    % Read TS
    CB  = load(fullfile(gray_dir,'Checkers','TSeries','Scan1','tSeries1.mat'));
    CBts = CB.tSeries;

    WE  = load(fullfile(gray_dir,'Words_English','TSeries','Scan1','tSeries1.mat'));
    WEts = WE.tSeries;

    WH  = load(fullfile(gray_dir,'Words_Hebrew','TSeries','Scan1','tSeries1.mat'));
    WHts = WH.tSeries;


    %% Filter to the ROI
    CBroi = CBts(:,ROIindices);
    CBecc_roi = CBecc(ROIindices);
    CBr2_roi = CBr2(ROIindices);
    CBsize_roi = CBsize(ROIindices);

    WEroi = WEts(:,ROIindices);
    WEecc_roi = WEecc(ROIindices);
    WEr2_roi = WEr2(ROIindices);
    WEsize_roi = WEsize(ROIindices);

    WHroi = WHts(:,ROIindices);
    WHecc_roi = WHecc(ROIindices);
    WHr2_roi = WHr2(ROIindices);
    WHsize_roi = WHsize(ROIindices);


    % FILTER FOR VE AND MIN ECC
    filterVE_ecc = find(CBr2_roi>VE & WEr2_roi>VE & WHr2_roi>VE ...
                    & CBecc_roi>min_ecc & WEecc_roi>min_ecc & WHecc_roi>min_ecc ...
                    & CBsize_roi>min_size & WEsize_roi>min_size & WHsize_roi>min_size);

    % Eliminate VEmin and fovea and too small size
    CBroi_VE_minecc = CBroi(:,filterVE_ecc);
    CBecc_roi_VE_minecc = CBecc_roi(filterVE_ecc);
    CBr2_roi_VE_minecc = CBr2_roi(filterVE_ecc);
    CBsize_roi_VE_minecc = CBsize_roi(filterVE_ecc);

    WEroi_VE_minecc = WEroi(:,filterVE_ecc);
    WEecc_roi_VE_minecc = WEecc_roi(filterVE_ecc);
    WEr2_roi_VE_minecc = WEr2_roi(filterVE_ecc);
    WEsize_roi_VE_minecc = WEsize_roi(filterVE_ecc);

    WHroi_VE_minecc = WHroi(:,filterVE_ecc);
    WHecc_roi_VE_minecc = WHecc_roi(filterVE_ecc);
    WHr2_roi_VE_minecc = WHr2_roi(filterVE_ecc);
    WHsize_roi_VE_minecc = WHsize_roi(filterVE_ecc);

    % Check
    filtered_voxels = ROIindices(filterVE_ecc);
    A = table((1:1:length(CBecc_roi_VE_minecc))', filtered_voxels, ...
        CBecc_roi_VE_minecc',WEecc_roi_VE_minecc',WHecc_roi_VE_minecc', ...
        CBsize_roi_VE_minecc',WEsize_roi_VE_minecc',WHsize_roi_VE_minecc');
    A.Properties.VariableNames = {'Index','Voxel','CBecc','WEecc','WHecc','CBsize','WEsize','WHsize'};
    A.Properties.VariableTypes = {'int32','int32','double','double','double','double','double','double'};

    % Normalize time series so that the values go from 0 to 1 using matlab's normalize, 'range' defaults to [0,1]
    CBnorm = normalize(CBroi_VE_minecc, 'range');
    WEnorm = normalize(WEroi_VE_minecc, 'range');
    WHnorm = normalize(WHroi_VE_minecc, 'range');

    % Obtain correlations for all timeseries
    CBWEcorr = diag(corr(CBnorm, WEnorm));
    CBWHcorr = diag(corr(CBnorm, WHnorm));
    WEWHcorr = diag(corr(WEnorm, WHnorm));

    if do_indiv_sub_plot

        % ind = filterVE_ecc(1);
        % ind = 3;
        indstr = input('Enter index voxel for time series: ', 's');
        ind = str2num(indstr);
        my_voxel = filtered_voxels(ind);

        titlestr1 = sprintf('ROI:%s | vox:%i | sub:%s | minVE:%i%% | minSize:%.2g | minEcc:%.2g', ...
            roi, my_voxel, strrep(subname,'_','\_'), round(100*VE), min_ecc, min_size);
        titlestr2 = sprintf(...
        '[ecc (R2) size] CB:%.2g (%i%%) %.2g | WE:%.2g (%i%%) %.2g | WH:%.2g (%i%%) %.2g', ...
        CBecc_roi_VE_minecc(ind), round(100*CBr2_roi_VE_minecc(ind)), CBsize_roi_VE_minecc(ind), ...
        WEecc_roi_VE_minecc(ind), round(100*WEr2_roi_VE_minecc(ind)), WEsize_roi_VE_minecc(ind), ...
        WHecc_roi_VE_minecc(ind), round(100*WHr2_roi_VE_minecc(ind)), WHsize_roi_VE_minecc(ind));

        % plot scatterplot
        hh = figure('visible','off');
        nrow = 4; ncol = 3;
        xmax = 10;; ymax = 10;
        subplot(nrow,ncol,1)
        plot(WEecc_roi_VE_minecc, CBecc_roi_VE_minecc, 'b.'); axis equal
        hold on; plot(WEecc_roi_VE_minecc(ind), CBecc_roi_VE_minecc(ind), 'ro')
        xlabel('ENGLISH'); xlim([0,xmax])
        ylabel('CHECKER'); ylim([0,ymax])
        identityLine(gca);
        title('Eccentricy [deg]')

        subplot(nrow,ncol,2)
        plot(WHecc_roi_VE_minecc, CBecc_roi_VE_minecc, 'b.'); axis equal
        hold on; plot(WHecc_roi_VE_minecc(ind), CBecc_roi_VE_minecc(ind), 'ro')
        xlabel('HEBREW'); xlim([0,xmax])
        ylabel('CHECKER'); ylim([0,ymax])
        identityLine(gca);

        subplot(nrow,ncol,3)
        plot(WEecc_roi_VE_minecc, WHecc_roi_VE_minecc, 'b.'); axis equal
        hold on; plot(WEecc_roi_VE_minecc(ind), WHecc_roi_VE_minecc(ind), 'ro')
        xlabel('ENGLISH'); xlim([0,xmax])
        ylabel('HEBREW'); ylim([0,ymax])
        identityLine(gca);

        subplot(nrow,ncol,4)
        hist(CBWEcorr);     
        xlabel('ENGLISH-CHECKERS');
        ylabel('COUNT');
        xlim([0,1]);
        title('Histogram of pairwise correlations')
        
        subplot(nrow,ncol,5)
        hist(CBWHcorr);
        xlim([0,1]);
        xlabel('HEBREW-CHECKERS');
    
        subplot(nrow,ncol,6)
        hist(WEWHcorr);
        xlim([0,1]);
        xlabel('ENGLISH-HEBREW');

        subplot(nrow,ncol,2*ncol+1:3*ncol)
        plot(2+CBnorm(:,ind), 'k-', 'LineWidth', 0.3); hold on;
        plot(1+WEnorm(:,ind), 'r-', 'LineWidth', 0.3);
        plot(WHnorm(:,ind), 'b-', 'LineWidth', 0.3);
        legend([{'CB'},{'WE'},{'WH'}]);
        % xlabel('Volumes'); 
        ylabel('Norm. BOLD [0,1]');
        grid on;
        title(titlestr1)

        subplot(nrow,ncol,3*ncol+1:4*ncol)
        % Obtain predictions
        CBmp = normalize(pmVistaObtainPrediction1voxel(CBresults, my_voxel), 'range');
        WEmp = normalize(pmVistaObtainPrediction1voxel(WEresults, my_voxel), 'range');
        WHmp = normalize(pmVistaObtainPrediction1voxel(WHresults, my_voxel), 'range');        
        plot(2+CBmp, 'k-', 'LineWidth', 0.3); hold on;
        plot(1+WEmp, 'r-', 'LineWidth', 0.3);
        plot(WHmp, 'b-', 'LineWidth', 0.3);
        legend([{'CB pred.'},{'WE pred.'},{'WH pred.'}])
        xlabel('Volumes'); 
        ylabel('Norm. BOLD [0,1]');
        grid on;
        title(titlestr2)

        fname = ['TS_roi-' roi '_vox-' num2str(ind) '-' num2str(my_voxel) '_sub-' subname ...
                 '_minVE-' num2str(VE) '_minSize-' num2str(min_ecc) '_minEcc-' num2str(min_size) ...
                 '_v2_resend-' res_end '.png'];
        fpath = fullfile(image_dir, fname)
        % Adjust margins
        set(hh, 'PaperUnits', 'inches');
        set(hh, 'PaperPosition', [0 0 10 10]); % Adjust the size as needed
        saveas(hh, fpath, 'png') 
    end

    if do_all_sub_plot
        all_CB_ecc    = [all_CB_ecc; CBecc_roi_VE_minecc'];
        all_WE_ecc    = [all_WE_ecc; WEecc_roi_VE_minecc'];
        all_WH_ecc    = [all_WH_ecc; WHecc_roi_VE_minecc'];
        all_CBWE_corr = [all_CBWE_corr; CBWEcorr];
        all_CBWH_corr = [all_CBWH_corr; CBWHcorr];
        all_WEWH_corr = [all_WEWH_corr; WEWHcorr];
    end
    
end

if do_all_sub_plot
    hh = figure('visible','off');
    nrow = 2; ncol = 3;
    xmax = 10; ymax = 10;

    subplot(nrow,ncol,1)
    plot(all_WE_ecc, all_CB_ecc, 'b.'); axis equal
    
    % c = ff_histogramHeat(all_WE_ecc, all_CB_ecc, [0,xmax], [0,ymax], 50);
    
                            %  cmapValuesHist,fov,roiName,fieldName,fontsize, cutoff)
    % ff_histogramHeat(x, y, maxValueX, maxValueY, numHistBins) 
    %
    % Makes a heat map!
    % INPUTS
    % x:            the vector of x values. assuming row vec!
    % y:            the vector of y values. assuming row vec!
    % maxValue:     limits of the x and y axis
    % numHistBins   something around 50
    % scatter(all_WE_ecc, all_CB_ecc, 'b.'); axis equal

    xlabel('ENGLISH'); xlim([0,xmax])
    ylabel('CHECKER'); ylim([0,ymax])
    identityLine(gca);
    title('Eccentricy [deg]')

    subplot(nrow,ncol,2)
    plot(all_WH_ecc, all_CB_ecc, 'b.'); axis equal
    xlabel('HEBREW'); xlim([0,xmax])
    ylabel('CHECKER'); ylim([0,ymax])
    identityLine(gca);

    subplot(nrow,ncol,3)
    plot(all_CB_ecc, all_WH_ecc, 'b.'); axis equal
    xlabel('ENGLISH'); xlim([0,xmax])
    ylabel('HEBREW'); ylim([0,ymax])
    identityLine(gca);

    subplot(nrow,ncol,4)
    hist(all_CBWE_corr);     
    xlabel('ENGLISH-CHECKERS');
    ylabel('COUNT');
    xlim([0,1]);
    title('Histogram of pairwise correlations')

    subplot(nrow,ncol,5)
    hist(all_CBWH_corr);
    xlim([0,1]);
    xlabel('HEBREW-CHECKERS');

    subplot(nrow,ncol,6)
    hist(WEWHcorr);
    xlim([0,1]);
    xlabel('ENGLISH-HEBREW');

    fname = ['GROUP_scatter_roi-' roi '_sub-' subname ...
             '_minVE-' num2str(VE) '_minSize-' num2str(min_ecc) ...
             '_minEcc-' num2str(min_size) ...
             '_v2_resend-' res_end '.png'];
    fpath = fullfile(image_dir, fname)
    % Adjust margins
    set(hh, 'PaperUnits', 'inches');
    set(hh, 'PaperPosition', [0 0 6 6]); % Adjust the size as needed
    saveas(hh, fpath, 'png') 
end
