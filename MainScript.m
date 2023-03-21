    %% (0) INIT
% This repository is based on the code created by Rosemary Le, part of coverageReading
% repository. For the Stimulus Dependence paper (2019) we tried to make it into
% a reproducible and reusable process. The thing is that we need to separate
% this specific project from others, we need to separate data from code.
tbUse PRF-StimDependence;

close all; clear all;
mrvCleanWorkspace;
cr         = struct();
cr.codeDir = sdRP;

% WHERE THE NEW DATA IS
cr.dirs.BASE     = '/black/localhome/glerma/TESTDATA/PRF-StimDependence';
cr.dirs.DATA     = fullfile(cr.dirs.BASE,'DATA');
cr.dirs.ANALYSIS = fullfile(cr.dirs.BASE,'ANALYSIS');
cr.dirs.ORG      = fullfile(cr.codeDir,'DATA','ANALYSIS','matlabfiles','organization');
cr.dirs.DEF      = fullfile(cr.codeDir,'DATA','ANALYSIS','matlabfiles','defineProjectDefaults');
% cr.dirs.FIG      = fullfile(cr.codeDir,'DATA','figures');
cr.dirs.FIG     = fullfile('/Users/glerma/Library/CloudStorage/GoogleDrive-garikoitz@gmail.com/My Drive/STANFORD/PROJECTS/2018 Reading across maps (Rosemary)/__PUBLISH__/2022_PNAS(3rd)','figures');
cr.dirs.FIGPNG  = fullfile(cr.dirs.FIG,'png');
cr.dirs.FIGSVG  = fullfile(cr.dirs.FIG,'svg');
% if ~isfolder(cr.dirs.FIG); mkdir(cr.dirs.FIG); end
% if ~isfolder(cr.dirs.FIGPNG); mkdir(cr.dirs.FIGPNG); end
% if ~isfolder(cr.dirs.FIGSVG); mkdir(cr.dirs.FIGSVG); end

% CONTINUE WITH THE NORMAL PROCESSING
% add to path the required matlab files inside the project, with info to run the project
% addpath(genpath(fullfile(cr.dirs.ANALYSIS,'matlabfiles')));

% Rosemary relied on this file that contains most of the subjects and other
% lists. Make it work with relative paths and store it in each project repository
% This file was used as well for:
% - copying files to a new location
% - editing mrSession to reflect the file changes
cr.bk = bookKeeping(cr);

%% (1) Run PRFs again
if 0
% subjects we want to do this for
list_subInds        = [31:36 38:44];  % Hebrew
% list_subInds      = [1:20];  % Original 20
% mw (13) for Words failed, continue with the next ones for now
% list_subInds      = [18:20];
%17 and 13 failed at beginning
% list_subInds     = [1,3,4,13:20];
% list_dtNames     = {'WordSmall','WordLarge'};
% list_dtNames     = {'Checkers'};

    for subind = list_subInds

        mrvCleanWorkspace;
        % subind  = list_subInds(ns);
        subname = cr.bk.list_sub{subind};
        [~,anatName]=fileparts(cr.bk.list_anatomy{subind});
        fprintf('\nSubDetails:\nInd:%i, StrInd:%s, subname:%s, Name:%s, anatName:%s\n',...
            subind,cr.bk.list_subNumberString{subind},subname,...
            cr.bk.list_names{subind},anatName)

        % Change dir, we need to run analysis where mrSession is
        % FOR ALL
        chdir(cr.bk.list_sessionRet{subind})
        prf.dirVistacc = cr.bk.list_sessionRet{subind};
        % FOR WORD LARGE SMALL
        % chdir(cr.bk.list_sessionSizeRet{subind})
        % prf.dirVistacc = cr.bk.list_sessionSizeRet{subind};

        %% PRF analysis
        % Read the generic params for all subjects
        run(fullfile(cr.dirs.DEF,'prfrun_defaults.m'));
        cr.defaults.prfrun.params = params;
        cr.defaults.prfrun.p      = p;
        clear('params'); clear('p');
        % Read prfRun_params specific to this subject
        % run(cr.bk.list_prfParams{subind}); NOT NECESSARY

        prf.dirAnatomy = cr.bk.list_anatomy{subind};
        prf.list_rmName= cr.bk.list_rmName{subind};
        prf.p.stimSize = cr.bk.list_stimSize(subind);
        prf.wSearch    = cr.bk.list_wSearch(subind);
        prf.prfModel   = cr.bk.list_prfModels{subind};
        prf.rois       = cr.bk.list_ROIs{subind};

        cr.subj.(subname).params.prf = prf;
        clear('prf');
        % This was on generics but requires specifics so... this is why I am
        % calling generics as many times as calling different subjects just in case
        cr.defaults.prfrun.params.stimSize = cr.subj.(subname).params.prf.p.stimSize;
        % Run the prfModel with mrVista
        % RUN USING mrVISTA NORMAL INSTALLATION
            cr = cr_prfRun(cr, subind);
            % Clean workspace of globals after each subject finishes
            mrvCleanWorkspace;
        % RUN USING container prfanalyze-vista:2.0.0 (no modelpred, we get r2)
            % Generate the config file
            % Run the container
            % pmLaunchDockerCommand('prfanalyze','ellipse','tr1dur300v3','afni6')
            % Convert the data back so that the rest of the scripts continue working
    end






end

%% -----------------------------------------------------------------------------
%% -----------------------------------------------------------------------------
%% (2) VAR EXPLAINED
inTheServer = false;
if inTheServer
    vw = initHiddenGray;
    vw = rmLoadDefault(vw);
    save('cc_vw.mat','vw')
else
    % Calculate variance explained between

    cd(fullfile(sdRP,'local'));
    % Load ROIs
    roiPath = fullfile(sdRP,'local','WangAtlas_V2v_left.mat');
    % roiPath = fullfile(sdRP,'local','lVOTRC.mat');
    AA      = load(roiPath);
    % vw   = initHiddenGray;
    % vw   = rmLoadDefault(vw);
    load(fullfile(sdRP,'local','cc_vw.mat'))
    vw    = loadROI(vw, roiPath, [],[],1,0);
    assert(isequal(AA.ROI.coords, vw.ROIs(1).coords));
    vw    = viewSet(vw, 'curdt', 'Checkers');
    V1indices   = viewGet(vw, 'roiIndices');
    % V1indices   = V1indices(1:1300);


    % Read checkers 1 and the average of 23
    C1  = load(fullfile(sdRP,'local','cc_checkers1_tSeries1.mat'));
    C23 = load(fullfile(sdRP,'local','cc_checkers23_tSeries1.mat'));
    C1  = C1.tSeries;
    C23 = C23.tSeries;

    % Filter to the ROI, left ventral V1
    C1V1      = C1(:,V1indices);
    C23V1     = C23(:,V1indices);
    subplot(3,3,1)
    meanC1V1  = mean(C1V1,2); plot(meanC1V1);hold on
    meanC23V1 = mean(C23V1,2); plot(meanC23V1,'r');
    title('Mean V1 activation');legend({'Checkers1','Checkers23'});
    % We need to demean it
    subplot(3,3,2)
    dmC1V1   = C1V1 - mean(C1V1);
    mdmC1V1  = mean(dmC1V1,2); plot(mdmC1V1);hold on
    dmC23V1  = C23V1 - mean(C23V1);
    mdmC23V1 = mean(dmC23V1,2); plot(mdmC23V1,'r');
    title('Demeaned mean V1 activation');legend({'Checkers1','Checkers23'});




    % Read the latest model fit
    C   = load(fullfile(sdRP,'local','retModel-Checkers-css-fFit-fFit.mat'));
    % Filter for V1
    rss    = C.model{1}.rss(V1indices);
    rawrss = C.model{1}.rawrss(V1indices);
    % Calculate var explained of model fit
    varexp         = 100*(1 - rss ./rawrss);
    subplot(3,3,4)
    plot(varexp,'.')
    ylim([0,100])
    title('var exp by index, new fit')
    subplot(3,3,5)
    histogram(varexp,100)
    xlim([0,100]);
    title('var exp by histogram, new fit')

    subplot(3,3,3)
    plot(dmC1V1,dmC23V1,'.');
    xlabel('dmC1V1')
    ylabel('dmC23V1')
    identityLine(gca)
    title('scatterplot of tSeries')

    subplot(3,3,6)
    histogram((diag(corr(dmC1V1,dmC23V1)).^2),200);
    xlabel('var explained')
    title('Corr of tSeries')



    % Calcualte var explained of the two time series
    F   = 100*(1 - sum((dmC1V1-dmC23V1).^2) ./ sum(dmC23V1.^2));
    F(F<0)=nan;
    subplot(3,3,7)
    plot(F,'.')
    ylim([0,100])
    title('var exp by index, 2 runs')
    subplot(3,3,8)
    histogram(F,100)
    xlim([0,100]);
    title('var exp by histogram, 2 runs')


    subplot(3,3,9)
    plot(varexp,F,'.'); xlim([0,100]); ylim([0,100]);identityLine(gca);
    xlabel('varexp fits')
    ylabel('varexp tSeries')
    title(sprintf('Medians: fit= %2.1f, tSeries= %2.1f', median(varexp,'omitnan'),median(F,'omitnan')))
end
%% -----------------------------------------------------------------------------

%% (3) PREPARE DATA: WORDS, CHECKERS AND FALSEFONTS
% Generate the rmroicell that we will use in all plots in this script
% This will read the results obtained by Rosemary or the re-run in 2021

readExisting = true;
whatFit = 'new';  % 'new' | 'Rosemary'

%
list_subInds  = [1:20];
% TEST what happens with dorsal
% V1-3d, V3a,IPS0-1
% list_roiNames = {'WangAtlas_V1v_left'
%                  'WangAtlas_V2v_left'
%                  'WangAtlas_V3v_left'
%                  'WangAtlas_hV4_left'
%                  'WangAtlas_VO1_left'
%                  'lVOTRC'
%                  'WangAtlas_IPS0'
%                  'WangAtlas_IPS1'};
list_roiNames = {'WangAtlas_V1d_left'
                 'WangAtlas_V2d_left'
                 'WangAtlas_V3d_left'
                 'WangAtlas_V3A_left'
                 'WangAtlas_IPS0'
                 'WangAtlas_IPS1'};
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

list_dtNames  = {'Checkers','Words','FalseFont'};
list_rmNames  = {'retModel-Checkers-css-fFit.mat'
                 'retModel-Words-css-fFit.mat'
                 'retModel-FalseFont-css-fFit.mat' };
%{
% Use the originals calculated by Rosemary
list_rmNames  = {'retModel-Checkers-css.mat'
                 'retModel-Words-css.mat'
                 'retModel-FalseFont-css.mat'};
%}
rmroiFname = ['rmroicell_subInds-1to20_dtNames-cb-w-ff_fits-' whatFit '_LeftRightROIs_2023.mat'];
rmroiFname='rmroicell_subInds-1to20_dtNames-cb-w-ff_fits-new_dorsalROIs_2023.mat';
if readExisting
    % Check if the file exists in the local dir, otherwise download from
    % OSF
    fpath = fullfile(sdRP,'local',rmroiFname);
    
    if ~isfile(fpath)
        % Download from OSF
        url = 'https://osf.io/download/y7wp6/';
        location = websave(fpath,url);   % >>> NOT WORKING
    end
    % Load it
    load(fpath,'rmroiCell');
    load(fullfile(sdRP,'DATA',rmroiFname))
    
    % If we have ventral ROIs, we should merge them with dorsal ones
    bu_rmroiCell = rmroiCell;
    bu_list_roiNames = list_roiNames;
    list_roiNames = {'WangAtlas_V1_left'
        'WangAtlas_V2_left'
        'WangAtlas_V3_left'
        'WangAtlas_hV4_left'
        'WangAtlas_VO1_left'
        'WangAtlas_V3A_left'
        'WangAtlas_IPS0_left'
        'WangAtlas_IPS1_left'
        'WangAtlas_V1_right'
        'WangAtlas_V2_right'
        'WangAtlas_V3_right'
        'WangAtlas_hV4_right'
        'WangAtlas_VO1_right'
        'WangAtlas_V3A_right'
        'WangAtlas_IPS0_right'
        'WangAtlas_IPS1_right'};
    rmroiCell = cell([size(bu_rmroiCell,1), (22-6), size(bu_rmroiCell,3)]);
    if strcmp('WangAtlas_V1v_left',bu_list_roiNames{4}) && length(bu_list_roiNames)==22
        nnr = 0;
        for nr=1:size(bu_rmroiCell,2)
            if ismember(nr,[1:3,12:14])
                
                nnr = nnr + 1
                for ii=1:3
                    for jj=1:size(bu_rmroiCell,1)
                        rmroiCell(jj,nnr,ii) = merge_cell_content(bu_rmroiCell,jj,nr,ii);
                    end
                end
            elseif ismember(nr,[7:11, 18:22])
                
                nnr = nnr + 1
                for ii=1:3
                    rmroiCell(:,nnr,ii) = bu_rmroiCell(:,nr,ii);
                end
            end
        end
        
        
        
    end
    
    
else
    rmroiCell = ff_rmroiCell(cr,...
        list_subInds,...
        list_roiNames,...
        list_dtNames, ...
        list_rmNames,...
        'list_path',cr.bk.list_sessionRet,...
        'latest_fFit',true, ...
        'checkYear','2022');
    % Save rmroicell
    save(fullfile(sdRP,'DATA',rmroiFname),'rmroiCell')
    
end

% Read the generic params for coverage for all subjects
cr.defaults.covfig.vfc = ff_vfcDefault();
cr.defaults.covfig.vfc.list_roiNames = list_roiNames;
% data types we want to look at
cr.defaults.covfig.vfc.list_dtNames = list_dtNames;
% names of the rm in each dt
cr.defaults.covfig.vfc.list_rmNames = list_rmNames;
% subinds = [31:36 38:44]; % Hebrew
% cr.defaults.covfig.vfc = ff_vfcDefault_Hebrew();

%% (4) Time series comparisons
readExisting = true;
whatFit  = 'new';
tsFname = ['tSeries_subInds-1to12_dtNames-w-ff_fits-' whatFit '.mat'];

if readExisting
    load(fullfile(sdRP,'DATA',tsFname),'tSs');
else

    rmroiCell_WFF    = rmroiCell(1:12,6,2:3);
    list_roiNames6   = list_roiNames(6);
    list_rmDescripts = {'Words','FalseFont'};
    tSs              = table();

    for subind=1:12
        subname = cr.bk.list_sub{subind};
        [~,anatName]=fileparts(cr.bk.list_anatomy{subind});
        fprintf('\nSubDetails:\nInd:%i, StrInd:%s, subname:%s, Name:%s, anatName:%s\n',...
            subind,cr.bk.list_subNumberString{subind},subname,...
            cr.bk.list_names{subind},anatName)
        % Select this subject
        thisW     = rmroiCell_WFF{subind,1,1};
        thisFF    = rmroiCell_WFF{subind,1,2};
        assert(isequal(thisW.indices, thisFF.indices))

        % Load time series
        Wts       = load(fullfile(cr.bk.list_sessionRet{subind},...
            'Gray','Words','TSeries','Scan1','tSeries1.mat'));
        Wts       = Wts.tSeries';
        FFts      = load(fullfile(cr.bk.list_sessionRet{subind},...
            'Gray','FalseFont','TSeries','Scan1','tSeries1.mat'));
        FFts      = FFts.tSeries';

        % Populate table
        tmpT      = table();
        tmpT.SUB  = categorical(repmat({subname},[length(thisW.indices),1]));
        tmpT.indx = thisW.indices;
        tmpT.W    = Wts(thisW.indices,:);
        tmpT.Wco  = 100*thisW.co';
        tmpT.Wecc = thisW.ecc';
        tmpT.Wx   = thisW.x0';
        tmpT.Wy   = thisW.y0';
        tmpT.Wsig = thisW.sigma';

        tmpT.FF   = FFts(thisFF.indices,:);
        tmpT.FFco = 100*thisFF.co';
        tmpT.FFecc= thisFF.ecc';
        tmpT.FFx  = thisFF.x0';
        tmpT.FFy  = thisFF.y0';
        tmpT.FFsig= thisFF.sigma';

        % Concatenate tables
        tSs = [tSs ; tmpT];
    end
    save(fullfile(sdRP,'DATA',tsFname),'tSs')
end


% Compare the timeSeries
% Filter: var expll > 20% & FFecc > 5deg
%
comin = 20;
tSsf = tSs(   tSs.FFco > comin & tSs.Wco > comin ...
            & tSs.FFecc > 5 & tSs.FFecc < 15 ...
            & tSs.Wecc  > 2 & tSs.Wecc  < 15 ...
            & tSs.SUB=='cc' &tSs.indx==1.1563e+05 ...
            , :);
            % & tSs.FFecc > tSs.Wecc ...

% Check by plotting it
xx = mrvNewGraphWin('check time series',[],true);
position = [0.005 0.062 .8 .8 ];
set(xx, 'position',position)
nrow = 2;
ncol = 2;


subplot(nrow,ncol,3)
plot(tSsf.Wecc, tSsf.FFecc,'ko');identityLine(gca)
xlabel('WORDS Eccentricity'); ylabel('FF Eccentricity');
set(gca,'FontSize',14)
axis equal
xlim([0,15]);ylim([0,15])

% Obtain matrices with the time series
W = tSsf.W;
F = tSsf.FF;

% Plot mean time series
subplot(nrow,ncol,1)
if size(W,1)>1
    plot(mean(W),'r-');hold on; plot(mean(F),'b-')
else
    plot(W,'r-');hold on; plot(F,'b-')
end
xlabel('time');legend({'Words','FalseFonts'}); ylabel('original signal')
title('Mean signals: original')
set(gca,'FontSize',14)

% Demean and plot again
dmW   = W - mean(W,2);
dmF   = F - mean(F,2);
subplot(nrow,ncol,2)
if size(dmW,1)>1
    plot(mean(dmW),'r-');hold on; plot(mean(dmF),'b-')
else
    plot(dmW,'r-');hold on; plot(dmF,'b-')
end
xlabel('time');legend({'Words','FalseFonts'});ylabel('demeaned signal')
title('Mean signals: demeaned')
set(gca,'FontSize',14)


% Substract and plot difference
if size(dmW,1)>1
    dmWF  = dmW - dmF;
    dmdiv = dmW ./ dmF;
    subplot(nrow,ncol,4)
    plot(mean(dmWF),'k-');hold on;plot(mean(dmdiv),'k--');
    xlabel('time');legend({'Mean Words-FalseFonts','Mean Words ./ FalseFonts'})
    title('Mean of demeaned W-FF'); ylabel('demeaned W-FF')
    set(gca,'FontSize',14)
else
    subplot(nrow,ncol,4)
    h1=viscircles([tSsf.Wx,tSsf.Wy],tSsf.Wsig,'Color','r','LineStyle','-'); hold on;
    h2=viscircles([tSsf.FFx,tSsf.FFy],tSsf.FFsig,'Color','b','LineStyle','-')
    plot([-15,15],[0,0],'k-','LineWidth',0.5)
    plot([0,0],[-15,15],'k-','LineWidth',0.5)
    legend([h1,h2],{'Words','FalseFonts'});

    xlabel('X'); ylabel('Y');
    title('pRF Center positions');
    set(gca,'FontSize',14)
    axis equal
    xlim([-15,15]);ylim([-15,15])
end


%% (5) Check median variance explained ver subject and ROI
% list_subInds  = [1:20];
% list_subInds  = [31:36 38:44];
subnames      = cr.bk.list_sub(list_subInds);
list_dtNames  = {'Checkers','Words'};
list_roiNames = {'WangAtlas_V1v_left'
                 'WangAtlas_V2v_left'
                 'WangAtlas_V3v_left'
                 'WangAtlas_hV4_left'
                 'WangAtlas_VO1_left'
                 'lVOTRC'};
rmroiCellMED  = rmroiCell(:,1:6,1:2);

% Create the table to store the results
howManyRows = 2 * length(list_subInds);
cols = {'sub','stim',list_roiNames{:},'nVoxVotrc'};
howManyCols = length(cols);
RES = array2table(nan(howManyRows,howManyCols));
% Change col names
RES.Properties.VariableNames = cols;
RES.sub  = categorical([subnames;subnames]);
RES.stim = categorical([repmat(list_dtNames(1),[length(list_subInds),1]);repmat(list_dtNames(2),[length(list_subInds),1])]);

for ii=1:length(list_subInds)
    nam = subnames{ii};
    for jj = 1:length(list_roiNames)
        roi  = list_roiNames{jj};
        for kk = 1:length(list_dtNames)
            dtName = list_dtNames{kk};
            % fprintf('\n\n%s >> %s >> %s\n',nam,roi,dtName)
            tmprss    = rmroiCell{ii,jj,kk}.rss;
            tmprawrss = rmroiCell{ii,jj,kk}.rawrss;
            R2        = 100*(1 - tmprss ./ tmprawrss);
            RES{RES.sub==nam & RES.stim==dtName,roi} = median(R2,'omitnan');
            if jj==6
                RES{RES.sub==nam & RES.stim==dtName,'nVoxVotrc'} = length(rmroiCell{ii,jj,kk}.indices);
            end
        end
    end
end

%% FIGURE 1: (C) Groups coverage plots

% (A) Explain how to obtain


% (B) Explain how to obtain


%
% With the new data the groups plots look different, but it seems that it
% is due to thresholds
% >> Check colormap limits  so that checker looks bigger than words

% Group COVERAGE plots, take all subjects from list_subInds

% Select subjects we want to plot
% subinds = [1:20]; % Stanford Subjects, 1 is gomez, find anatomicals
% subinds = [1:12,14:16,18:20];

%{
rmroiCell_WC     = rmroiCell(:,1:6,1:2);
rmroiCell_WC     = flip(rmroiCell_WC,3);
list_roiNames16  = list_roiNames(1:6);
%}
% {
rmroiCell_VOTRC     = rmroiCell(:,6,1:2);
rmroiCell_VOTRC     = flip(rmroiCell_VOTRC,3);
list_roiNamesVOTRC  = list_roiNames(6);
list_dtNamesWC      = list_dtNames([2,1]);
list_rmNamesWC      = list_rmNames([2,1]);
%}

% Launch the function
fname = 'CoverageBoot_';%'Fig1_'; % '' for not saving
% fname = ''; %
[RF_mean, RF_individuals,empties] = figFunction_coverage_maxProfile_group(cr,list_subInds, ...
                                      'flip',false, ...
                                      'bootcontour', false, ...
                                      'rmroiCell',rmroiCell_VOTRC, ...
                                      'list_roiNames', list_roiNamesVOTRC, ...
                                      'list_dtNames', list_dtNamesWC, ...
                                      'list_rmNames', list_rmNamesWC, ...
                                      'sizedegs',15,...
                                      'minvarexp', 0.2, ...
                                      'numboots',25, ...
                                      'fname', fname, ...
                                      'vers',['v02_' whatFit 'fit'],...
                                      'invisible',false);

% PLOT THEM FOR VOTRC, DO BOOTSTRAPPING AND AVERAGE IT
Wind = [1:17,19:20];
ALLW=RF_individuals{1}(:,:,Wind);
ALLC=RF_individuals{2};

MVALS = zeros(128,128,50);
SVALS = zeros(128,128,50);
DVALS = zeros(128,128,50);

for kk=1:50
    % Remove 1 each time and create same plots with the remaining one
    randReplacement = datasample(1:9,9);
    alleng          = ALLW(:,:,randReplacement);
    allheb          = ALLC(:,:,randReplacement);

    % Calculate measures
    mval   = mean(allheb - alleng, 3);
    stdval = std(allheb - alleng, [],3);
    Cd     = zeros(128,128);
    for ii=1:128;for jj=1:128
            Cd(ii,jj)=computeCohen_d(allheb(ii,jj,:),alleng(ii,jj,:),'paired');
    end;end

    % Accummlate it
    MVALS(:,:,kk) = mval;
    SVALS(:,:,kk) = stdval;
    DVALS(:,:,kk) = Cd;

end

% obtain means again
mval   = mean(MVALS,3);
stdval = mean(SVALS,3);
Cd     = mean(DVALS,3);


% PLOT
mrvNewGraphWin('bootstrapsEngvsHebFOV','wide',true);

% subplot(1,3,1)
% % subplot(1,2,1)
% imagesc(mval);axis equal;colormap(parula);colorbar;grid
% title('MEAN: mean of 50 bootstraps [CB-Word]')
% xlim([1,128]);ylim([1,128])
% xticks([1,64,128]); yticks([1,64,128])
% xticklabels([-15,0,15]); yticklabels([-15,0,15])
% xlabel('Degs'); ylabel('Degs')
% % subplot(1,2,2)
% % surf(mval)


% subplot(1,3,2)
% imagesc(stdval);axis equal;colormap(parula);colorbar;grid
% title('SD: mean 50 crossvals [CB-Word]')
% xlim([1,128]);ylim([1,128])
% xticks([1,64,128]); yticks([1,64,128])
% xticklabels([-15,0,15]); yticklabels([-15,0,15])
% xlabel('Degs'); ylabel('Degs')

% subplot(1,3,3)
subplot(1,2,1)
imagesc(Cd);axis equal;colormap(parula);colorbar;grid
title("d': mean 50 bootstraps [CB-Word]")
xlim([1,128]);ylim([1,128])
xticks([1,64,128]); yticks([1,64,128])
caxis([-1.25,1.25])
xticklabels([-15,0,15]); yticklabels([-15,0,15])
xlabel('Degs'); ylabel('Degs')

% subplot(1,2,2)
subplot(1,2,2)
[X,Y] = meshgrid(1:128,1:128);
XX = ((X-64)/64)*15;
YY = ((Y-64)/64)*15;
YY = flipud(YY);
surf(XX,YY,Cd);
xlabel('X (degs)'); ylabel('Y (degs)')
zlabel("Cohen's d")
xlim([-15,15])
ylim([-15,15])
zlim([-0.8,0.8])
xticks([-15,-10,-5,0,5,10,15])
yticks([-15,-10,-5,0,5,10,15])
xticklabels({'-15','-10','-5','0','5','10','15'})
yticklabels({''})
set(gca,'FontSize',18)

fname = 'FOV_Comparisons_WordsvsCB';
set(0, 'DefaultFigureRenderer', 'painters');
saveas(gcf, fullfile(cr.dirs.FIGPNG, [fname '.png']), 'png')
saveas(gcf, fullfile(cr.dirs.FIGSVG,[fname '.svg']), 'svg')

%% FIGURE 2: (C) Eccentricity: Scatterplots: word-checkerboard (no IPS)
% Order is CB, W, FF, invert it so that it is W then CB

% LEFT

rmroiCell_WC     = rmroiCell(:,1:6,1:3);
rmroiCell_WC     = flip(rmroiCell_WC,3);
list_roiNames16  = list_roiNames(1:6);
list_rmDescripts = {'FalseFont', 'Checkers'};%     {'FalseFont'}


% Obtain equally thresholded voxels to scatterplot
varExplained=0.2;
[R,C_data,cr]=crThreshGetSameVoxel(cr,...
                                   rmroiCell_WC,...
                                   list_subInds,...
                                   list_roiNames16,...
                                   'cothres', varExplained,...
                                   'fieldrange', 15);

% Plot it
fontsize = 12;
fname = ['LEFT_scatterplot_eccentricity_FFVsCheck_6DorsalROIs_20subs_' whatFit 'Fit_v01'];
[percAboveIdentity,~] = crCreateScatterplot(R,C_data,cr,...
                                    list_subInds,...
                                    list_roiNames16,...
                                    list_rmDescripts,...
                                    'ecc', ...  % 'co'
                                    fontsize, ...
                                    '');

% RIGHT

rmroiCell_WC     = rmroiCell(:,9:16,1:2);
rmroiCell_WC     = flip(rmroiCell_WC,3);
list_roiNames16  = list_roiNames(9:16);
list_rmDescripts = {'Words','Checkers'};%    {'FalseFont'}


% Obtain equally thresholded voxels to scatterplot
varExplained=0.2;
[R,C_data,cr]=crThreshGetSameVoxel(cr,...
                                   rmroiCell_WC,...
                                   list_subInds,...
                                   list_roiNames16,...
                                   'cothres', varExplained,...
                                   'fieldrange', 15);

% Plot it
fontsize = 12;
fname = ['RIGHT_scatterplot_eccentricity_WordVsCheck_6DorsalROIs_20subs_' whatFit 'Fit_v01'];
[percAboveIdentity,~] = crCreateScatterplot(R,C_data,cr,...
                                    list_subInds,...
                                    list_roiNames16,...
                                    list_rmDescripts,...
                                    'ecc', ...  % 'co'
                                    fontsize, ...
                                    '');

                                
                                
                                
                                
%% FIGURE 3: (B) Line plots: word-checkerboard
% Uses the same voxel calculations from the previous plot
% If doubt or this is moved, calculate it again here

% Plot it
A = colormap(jet);close all;
fname = ['lineplot_WordVsCheck_6ROIs_20subs_' whatFit 'Fit_v02'];
% fname = '';
[diffs15,diffs5,posangles15,negangles15,posangles5,negangles5] = ...
                         crCreateLinePlot(R,cr,list_roiNames16,fname);





% plot angles
A = A(1:round(length(A)/length(list_roiNames16)):end,:);
Y = cell(2,2);
Y(1,:) = {posangles5,negangles5};
Y(2,:) = {posangles15,negangles15};

xx = mrvNewGraphWin('posangles15','wide');
set(xx,'Position',[0.005 0.062 .1 .3]);
for nx=1:2
    subplot(2,1,nx)
    Hs = [];
    legs = {};
    for ii=1:length(list_roiNames16)
        % subplot(2,3,ii)
        X = [Y{nx,1}{ii},Y{nx,2}{ii}];
        % histogram(diffs{ii}); hold on
        if ii==6
            H = histfit(X,20,'kernel'); hold on
        else
            H = histfit(X,100,'kernel'); hold on
        end
        H(1).FaceColor=A(ii,:);%[1,1,1];
        H(1).EdgeColor=A(ii,:);%[1,1,1];
        H(1).FaceAlpha=0;
        H(1).EdgeAlpha=0;
        H(2).Color=A(ii,:);
        H(2).LineWidth=1.5;
        H(2).Visible='on';
        H(1).YData=H(1).YData/(sum(H(1).YData) * (H(1).XData(2) - H(1).XData(1)) );
        H(2).YData=H(2).YData/trapz(H(2).XData, H(2).YData);
        Hs = [Hs, H(2)];
        % plot([0,0],[0,0.1],'k-.')
        % ylim([0,0.1])
        % xlim([-90,90])
        % disp(trapz(H(2).XData, H(2).YData))
        % leg  = strrep(strrep(list_roiNames16{ii},'WangAtlas_',''),'_','\_');
        leg  = strrep(strrep(list_roiNames16{ii},'WangAtlas_',''),'_left','');
        leg  = strrep(leg,'l','');
        leg  = [leg ' (N=' num2str(length(X)) ')'];
        legs = [legs {leg}];
    end
    ylim([0,0.04])
    xlim([-90,90])
    plot([0,0],[0,0.1],'k-.')
    set(gca,'FontSize', 8)
    legend(Hs,legs,'FontSize', 8)
    if nx==1; title('Perifoveal T','FontSize', 10);end
    if nx==2; title('Near periphery T','FontSize', 10);end
    xlabel('Angle [deg]','FontSize', 10)
    ylabel('Relative count','FontSize', 10)

end

% SAVE THE FIG
set(gcf,'color','w');
fname = ['lineplotANGLES_WordVsCheck_6ROIs_20subs_' whatFit 'Fit_v02'];
if ~isempty(fname)
    saveas(gcf, fullfile(cr.dirs.FIGPNG, [fname '.png']), 'png')
    saveas(gcf, fullfile(cr.dirs.FIGSVG,[fname '.svg']), 'svg')
end










% Plot distances
Y = cell(2,1);
Y(1,1) = {diffs5};
Y(2,1) = {diffs15};

xx = mrvNewGraphWin('Differences','wide');
set(xx,'Position',[0.005 0.062 .1 .2 ]);
for nx=1:2
    subplot(2,1,nx)
    Hs = [];
    legs = {};
    for ii=1:length(list_roiNames16)
        % subplot(2,3,ii)
        X = [Y{nx,1}{ii}];
        H = histfit(X,100,'kernel'); hold on
        H(1).FaceColor=A(ii,:);
        H(1).EdgeColor=A(ii,:);
        H(1).FaceAlpha=0;
        H(1).EdgeAlpha=0;
        H(2).Color=A(ii,:);
        H(2).LineWidth=1.5;
        H(2).Visible='on';
        H(1).YData=H(1).YData/(sum(H(1).YData) * (H(1).XData(2) - H(1).XData(1)) );
        H(2).YData=H(2).YData/trapz(H(2).XData, H(2).YData);
        Hs = [Hs, H(2)];
        % plot([0,0],[0,0.5],'k-.')
        % ylim([0,0.1])
        % xlim([-5,10])
        leg  = strrep(strrep(list_roiNames16{ii},'WangAtlas_',''),'_left','');
        leg  = strrep(leg,'l','');
        leg  = [leg ' (N=' num2str(length(X)) ')'];
        legs = [legs {leg}];
    end
    ylim([0,1])
    xlim([-4,8])
    % if nx==1;xlim([-4,4]);end
    % if nx==2;xlim([-5,10]);end
    plot([0,0],[0,1],'k-.')
    set(gca,'FontSize', 8)
    legend(Hs,legs,'FontSize', 8)
    if nx==1; title('Perifoveal D','FontSize', 10);end
    if nx==2; title('Near periphery D','FontSize', 10);end
    xlabel('Delta Eccentricity Checkers - Words [deg]')
    ylabel('Relative count')
end

% SAVE THE FIG
set(gcf,'color','w');
fname = ['lineplotDIFF_WordVsCheck_6ROIs_20subs_' whatFit 'Fit_v02'];
if ~isempty(fname)
    saveas(gcf, fullfile(cr.dirs.FIGPNG, [fname '.png']), 'png')
    saveas(gcf, fullfile(cr.dirs.FIGSVG,[fname '.svg']), 'svg')
end

%

%% FIGURE 4: (A) Variance Explained: Scatterplot: word-checkerboard
% Uses the same voxel calculations from the previous plot
% If doubt or this is moved, calculate it again here

% Plot it
fname=['scatterplot_varianceExplained_WordVsCheck_6ROIs_20subs_' whatFit '_Fit_v02'];
fname='';
[percAboveIdentity,perROI] = crCreateScatterplot(R,C_data,cr,...
                             list_subInds,...
                             list_roiNames16,...
                             list_rmDescripts,...
                             'co', ...  % 'co'
                             fontsize, ...
                             fname);


%% FIGURE 5: (A) Eccentricity and (B) Variance Explained: Scatterplot: word-falsefont
% Order is CB, W, FF, invert it so that it is W then CB
rmroiCell_WF    = rmroiCell(:,1:6,2:3);
list_roiNames16 = list_roiNames(1:6);
list_rmDescripts = {'Words','FalseFont'};
fname='';  % for not creating output
% Obtain equally thresholded voxels to scatterplot
[R,C_data,cr]=crThreshGetSameVoxel(cr,...
                                   rmroiCell_WF,...
                                   list_subInds,...
                                   list_roiNames16,...
                                   'cothres', 0.2,...
                                   'fieldrange', 15);
% FIG3A and FIGS2a
fname = ['scatterplot_eccentricity_WordVsFF_6ROIs_20subs_' whatFit 'Fit_v02'];
fname='';
crCreateScatterplot(R,C_data,cr,...
                    list_subInds,...
                    list_roiNames16,...
                    list_rmDescripts,...
                    'ecc', ...  % 'co'
                    fontsize,...
                    fname);


% % FIGS2b: Variance Explained: Scatterplot: word-falsefont
fname = ['scatterplot_varianceExplained_WordVsFF_6ROIs_20subs_' whatFit 'Fit_v02'];
fname='';
crCreateScatterplot(R,C_data,cr,...
                    list_subInds,...
                    list_roiNames16,...
                    list_rmDescripts,...
                    'co', ...  % 'co'
                    fontsize,...
                    fname);


% % FIG3B and FIGS3a: Eccentricity: Scatterplot: FF-checkerboard
rmroiCell_FC    = rmroiCell(:,1:6,[1,3]);
rmroiCell_FC    = flip(rmroiCell_FC,3);
list_roiNames16 = list_roiNames(1:6);
list_rmDescripts = {'FalseFont','Checkers'};

% Obtain equally thresholded voxels to scatterplot
[R,C_data,cr]=crThreshGetSameVoxel(cr,...
                                   rmroiCell_FC,...
                                   list_subInds,...
                                   list_roiNames16,...
                                   'cothres', 0.2,...
                                   'fieldrange', 15);

fname = ['scatterplot_eccentricity_FFVsCB_6ROIs_20subs_' whatFit 'Fit_v02'];
fname = '';
crCreateScatterplot(R,C_data,cr,...
                    list_subInds,...
                    list_roiNames16,...
                    list_rmDescripts,...
                    'ecc', ...  % 'co'
                    fontsize,...
                    fname);
% FIG S3b
fname = ['scatterplot_varianceExplained_FFVsCB_6ROIs_20subs_' whatFit 'Fit_v02'];
crCreateScatterplot(R,C_data,cr,...
                    list_subInds,...
                    list_roiNames16,...
                    list_rmDescripts,...
                    'co', ...  % 'co'
                    fontsize,...
                    fname);

%% FIGURE S1: ALL Eccentricity: Scatterplots: word-checkerboard
% Order is CB, W, FF, invert it so that it is W then CB
% (I run (0) and (3) before this if starting here
rmroiCell_WCIPS     = rmroiCell(:,:,1:2);
rmroiCell_WCIPS     = flip(rmroiCell_WCIPS,3);
list_rmDescripts    = {'Words','Checkers'};

ves=[0.20,0.05];
for ve=ves
    % Obtain equally thresholded voxels to scatterplot
    [R,C_data,cr]=crThreshGetSameVoxel(cr,...
                                       rmroiCell_WCIPS,...
                                       list_subInds,...
                                       list_roiNames,...
                                       'cothres', ve,...
                                       'fieldrange', 15);

    % Plot it
    fontsize = 20;
    fname = ['scatterplot_eccentricity_WordVsCB_IPS_VE' num2str(100*ve) '_20subs_' whatFit 'Fit_v02'];
    crCreateScatterplot(R,C_data,cr,...
                        list_subInds,...
                        list_roiNames,...
                        list_rmDescripts,...
                        'ecc', ...
                        fontsize ,...
                        fname);
    fname = ['scatterplot_varianceExplained_WordVsCheck_IPS_VE' num2str(100*ve) '_20subs_' whatFit 'Fit_v02'];
    crCreateScatterplot(R,C_data,cr,...
                        list_subInds,...
                        list_roiNames,...
                        list_rmDescripts,...
                        'co', ...
                        fontsize,...
                        fname);
end

%% FIGURE S2: ALL Eccentricity: Scatterplots: word-false font
% Order is CB, W, FF, invert it so that it is W then CB
rmroiCell_WCIPS     = rmroiCell(:,:,[2,3]);
list_rmDescripts    = {'Words','FalseFonts'};

ves=[0.20,0.05];
for ve=ves
    % Obtain equally thresholded voxels to scatterplot
    [R,C_data,cr]=crThreshGetSameVoxel(cr,...
                                       rmroiCell_WCIPS,...
                                       list_subInds,...
                                       list_roiNames,...
                                       'cothres', ve,...
                                       'fieldrange', 15);

    % Plot it
    fontsize = 20;
    fname = ['scatterplot_eccentricity_WordVsFF_IPS_VE' num2str(100*ve) '_20subs_' whatFit 'Fit_v02'];
    crCreateScatterplot(R,C_data,cr,...
                        list_subInds,...
                        list_roiNames,...
                        list_rmDescripts,...
                        'ecc', ...
                        fontsize ,...
                        fname);
    fname = ['scatterplot_varianceExplained_WordVsFF_IPS_VE' num2str(100*ve) '_20subs_' whatFit 'Fit_v02'];
    crCreateScatterplot(R,C_data,cr,...
                        list_subInds,...
                        list_roiNames,...
                        list_rmDescripts,...
                        'co', ...
                        fontsize,...
                        fname);
end

%% FIGURE S3: ALL Eccentricity: Scatterplots: FF-checkerboard
% Order is CB, W, FF, invert it so that it is W then CB
rmroiCell_WCIPS     = rmroiCell(:,:,[1,3]);
rmroiCell_WCIPS     = flip(rmroiCell_WCIPS,3);
list_rmDescripts    = {'FalseFonts','Checkers'};

ves=[0.20,0.05];
for ve=ves
    % Obtain equally thresholded voxels to scatterplot
    [R,C_data,cr]=crThreshGetSameVoxel(cr,...
                                       rmroiCell_WCIPS,...
                                       list_subInds,...
                                       list_roiNames,...
                                       'cothres', ve,...
                                       'fieldrange', 15);

    % Plot it
    fontsize = 20;
    fname = ['scatterplot_eccentricity_CBVsFF_IPS_VE' num2str(100*ve) '_20subs_' whatFit 'Fit_v02'];
    crCreateScatterplot(R,C_data,cr,...
                        list_subInds,...
                        list_roiNames,...
                        list_rmDescripts,...
                        'ecc', ...
                        fontsize ,...
                        fname);
    fname = ['scatterplot_varianceExplained_CBVsFF_IPS_VE' num2str(100*ve) '_20subs_' whatFit 'Fit_v02'];
    crCreateScatterplot(R,C_data,cr,...
                        list_subInds,...
                        list_roiNames,...
                        list_rmDescripts,...
                        'co', ...
                        fontsize,...
                        fname);
end

%% -----------------------------------------------------------------------------

%% PREPARE DATA: WORDS LARGE AND SMALL
% Generate the rmroicell that we will use in all plots in this script
% This will read the results obtained by Rosemary or the re-run in 2021

readExisting = true;
% Do the same with the small and large words
list_subInds     = [1,3,4,13:20];
list_dtNames     = {'WordSmall','WordLarge'};
list_roiNames = {'WangAtlas_V1v_left'
                 'WangAtlas_V2v_left'
                 'WangAtlas_V3v_left'
                 'WangAtlas_hV4_left'
                 'WangAtlas_VO1_left'
                 'lVOTRC'
                 'WangAtlas_IPS0'
                 'WangAtlas_IPS1'};
list_rmNames     = {'retModel-WordSmall-css.mat'
                    'retModel-WordLarge-css.mat' };
list_rmDescripts = {'WordSmall'...
                    'WordLarge'};
if readExisting
    load(fullfile(sdRP,'DATA',...
      'rmroicell_subInds-1-3-4-13to20_dtNames-Wsmall-Wlarge_fits-Rosemary.mat'),'rmroiCell');
else
    rmroiCell=ff_rmroiCell(cr,list_subInds,list_roiNames,list_dtNames,...
                           list_rmNames,'list_path',cr.bk.list_sessionSizeRet);
    % Save rmroicell just in case
    save(fullfile(sdRP,'DATA',...
      'rmroicell_subInds-1-3-4-13to20_dtNames-Wsmall-Wlarge_fits-Rosemary.mat'),'rmroiCell')
end

% Read the generic params for coverage for all subjects
cr.defaults.covfig.vfc = ff_vfcDefault();
cr.defaults.covfig.vfc.list_roiNames = list_roiNames;
% data types we want to look at
cr.defaults.covfig.vfc.list_dtNames = list_dtNames;
% names of the rm in each dt
cr.defaults.covfig.vfc.list_rmNames = list_rmNames;
% subinds = [31:36 38:44]; % Hebrew
% cr.defaults.covfig.vfc = ff_vfcDefault_Hebrew();

%% FIGURE xxx: (A) Eccentricity and (B) Variance Explained: Scatterplot: 'WordSmall','WordLarge'
% Order is WS,WL
list_roiNames16 = list_roiNames(1:6);
list_rmDescripts = {'WordSmall','WordLarge'};

% Obtain equally thresholded voxels to scatterplot
[R,C_data,cr]=crThreshGetSameVoxel(cr,...
                                   rmroiCell,...
                                   list_subInds,...
                                   list_roiNames16,...
                                   'cothres', 0.2,...
                                   'fieldrange', 15);

fname = 'scatterplot_eccentricity_WordSmallVsWordLarge_6ROIs_11subs_RosemaryFit_v01';
crCreateScatterplot(R,C_data,cr,...
                    list_subInds,...
                    list_roiNames16,...
                    list_rmDescripts,...
                    'ecc', ...  % 'co'
                    fname);


% FIGURE 6: (B) Variance Explained: Scatterplot: wordLarge-WordSmall
% Uses the same voxel calculations from the previous plot
% If doubt or this is moved, calculate it again here

% Plot it
fname = 'scatterplot_varianceExplained_WordSmallVsWordLarge_6ROIs_11subs_RosemaryFit_v01';
crCreateScatterplot(R,C_data,cr,...
                    list_subInds,...
                    list_roiNames16,...
                    list_rmDescripts,...
                    'co', ...  % 'co'
                    fname);

%% -----------------------------------------------------------------------------

%% (6) PREPARE DATA: ENGLISH AND HEBREW WORDS
% Generate the rmroicell that we will use in all plots in this script
% This will read the results obtained by Rosemary or the re-run in 2021
readExisting = true;
% Do the same with the small and large words
list_subInds     = [31:36 38:44];
list_roiNames = {'WangAtlas_V1v_left'
                 'WangAtlas_V2v_left'
                 'WangAtlas_V3v_left'
                 'WangAtlas_hV4_left'
                 'WangAtlas_VO1_left'
                 'lVOTRC'
                 'WangAtlas_IPS0'
                 'WangAtlas_IPS1'};

whatFit = 'new';  % 'new' | 'Rosemary'
list_dtNames     = {'Words_English','Words_Hebrew'};
list_rmDescripts = {'Words_English','Words_Hebrew'};
if strcmp(whatFit,'Rosemary')
    list_rmNames     = {'retModel-Words_English-css.mat','retModel-Words_Hebrew-css.mat' };
else
     list_rmNames     = {'retModel-Words_English-css-fFit.mat','retModel-Words_Hebrew-css-fFit.mat' };
end
matname = ['rmroicell_subInds-31to36-38to44_dtNames-WE-WH_fits-' whatFit '.mat'];
matname = ['rmroicell_subInds-31to36-38to44_dtNames-ALL_fits-' whatFit '_2022.mat'];

if readExisting
    load(fullfile(sdRP,'DATA',matname),'rmroiCell');
else
    rmroiCell=ff_rmroiCell(cr,list_subInds,list_roiNames,list_dtNames,...
                           list_rmNames,'list_path',cr.bk.list_sessionRet,...
                             'latest_fFit',true, ...
                             'checkYear','2022');
    % Save rmroicell just in case
    save(fullfile(sdRP,'DATA',matname),'rmroiCell')
end

% Read the generic params for coverage for all subjects
cr.defaults.covfig.vfc = ff_vfcDefault();
cr.defaults.covfig.vfc.list_roiNames    = list_roiNames;
% data types we want to look at
cr.defaults.covfig.vfc.list_dtNames     = list_dtNames;
% names of the rm in each dt
cr.defaults.covfig.vfc.list_rmNames     = list_rmNames;
cr.defaults.covfig.vfc.list_rmDescripts = list_rmDescripts;
% subinds = [31:36 38:44]; % Hebrew
% cr.defaults.covfig.vfc = ff_vfcDefault_Hebrew();

%% FIGURE 4: (A) Eccentricity and (B) Variance Explained: Scatterplot: EW-HW
% Order is WS,WL
list_roiNames16 = list_roiNames(1:6);
list_rmDescripts = {'Words_English','Words_Hebrew'};

% Obtain equally thresholded voxels to scatterplot
cothresh = 0.2;
[R,C_data,cr]=crThreshGetSameVoxel(cr,...
                                   rmroiCell,...
                                   list_subInds,...
                                   list_roiNames16,...
                                   'cothres', cothresh,...
                                   'fieldrange', 7);


fprintf('------------------------------------\n')
fprintf('         %s - %s (R2:%g)        \n', list_rmDescripts{2},list_rmDescripts{1},cothresh)
fprintf('------------------------------------\n\n')
for ii=1:6
    fprintf('%s(N=%i)\n',strrep(list_roiNames16{ii},'WangAtlas_',''),length(R.X_rm2{ii}))
    [H P CI] = ttest(R.X_rm2{ii}-R.X_rm1{ii});
    fprintf('(X) P:%g, CI: [%g %g]\n',P,CI(1),CI(2))

    [H P CI] = ttest(R.Y_rm2{ii}-R.Y_rm1{ii});
    fprintf('(Y) P:%g, CI: [%g %g]\n\n',P,CI(1),CI(2))
end


fname = ['scatterplot_eccentricity_WordEngVsWordHeb_6ROIs_13subs_' whatFit 'Fit_v02'];
% fname = '';
crCreateScatterplot(R,C_data,cr,...
                    list_subInds,...
                    list_roiNames16,...
                    list_rmDescripts,...
                    'ecc', ...  % 'co'
                    fontsize, ...
                    fname);

crCreateScatterplot(R,C_data,cr,...
                    list_subInds,...
                    list_roiNames16,...
                    list_rmDescripts,...
                    'x0', ...  % 'co'
                    fontsize, ...
                    fname);

crCreateScatterplot(R,C_data,cr,...
                    list_subInds,...
                    list_roiNames16,...
                    list_rmDescripts,...
                    'y0', ...  % 'co'
                    fontsize, ...
                    fname);

% FIGURE 6: (B) Variance Explained: Scatterplot: wordLarge-WordSmall
% Uses the same voxel calculations from the previous plot
% If doubt or this is moved, calculate it again here

% Plot it
fname = ['scatterplot_varianceExplained_WordEngVsWordHeb_6ROIs_13subs_' whatFit 'Fit_v02'];
% fname = '';
crCreateScatterplot(R,C_data,cr,...
                    list_subInds,...
                    list_roiNames16,...
                    list_rmDescripts,...
                    'co', ...  % 'co', 'ecc'
                    fontsize, ...
                    fname);


%% FIGURE S5: Eccentricity and VE: Scatterplots: word-Eng and Heb

list_rmDescripts = {'Words_English','Words_Hebrew'};

ves=[0.20,0.05];
for ve=ves
    % Obtain equally thresholded voxels to scatterplot
    [R,C_data,cr]=crThreshGetSameVoxel(cr,...
                                       rmroiCell(:,:,[1,2]),...
                                       list_subInds,...
                                       list_roiNames,...
                                       'cothres', ve,...
                                       'fieldrange', 7);

    % Plot it
    fontsize = 20;% fname='';
    fname = ['scatterplot_eccentricity_WEngVsWHeb_IPS_VE' num2str(100*ve) '_13subs_' whatFit 'Fit_v02'];
    crCreateScatterplot(R,C_data,cr,...
                        list_subInds,...
                        list_roiNames,...
                        strrep(list_rmDescripts,'_','\_'),...
                        'ecc', ...
                        fontsize ,...
                        fname);
    fname = ['scatterplot_varianceExplained_WordVsCheck_IPS_VE' num2str(100*ve) '_13subs_' whatFit 'Fit_v02'];
    crCreateScatterplot(R,C_data,cr,...
                        list_subInds,...
                        list_roiNames,...
                        list_rmDescripts,...
                        'co', ...
                        fontsize,...
                        fname);
end

 %% FIGURE S6: X and Y: Scatterplots: word-Eng and Heb

list_rmDescripts = {'Words_English','Words_Hebrew'};

ve=0.20;
meass={'x0','y0'};
for meas=meass
    % Obtain equally thresholded voxels to scatterplot
    [R,C_data,cr]=crThreshGetSameVoxel(cr,...
                                       rmroiCell(:,:,[1,2]),...
                                       list_subInds,...
                                       list_roiNames,...
                                       'cothres', ve,...
                                       'fieldrange', 7);

    % Plot it
    fontsize = 20;% fname='';
    fname = ['scatterplot_' meas{:} '_WEngVsWHeb_IPS_VE' num2str(100*ve) '_13subs_' whatFit 'Fit_v02'];
    crCreateScatterplot(R,C_data,cr,...
                        list_subInds,...
                        list_roiNames,...
                        strrep(list_rmDescripts,'_','\_'),...
                        meas{:}, ...
                        fontsize ,...
                        fname);

end



%% FIGURE S7: Groups coverage plots: EW-HW
% RUN 0 and 6 before this


%{
% With the new data the groups plots look different, but it seems that it
% is due to thresholds
% >> Check colormap limits  so that checker looks bigger than words

% Group COVERAGE plots, take all subjects from list_subInds
%{
rmroiCell_noIPS = rmroiCell(:,1:6,1:2);
list_roiNames   = {'WangAtlas_V1v_left'
                   'WangAtlas_V2v_left'
                   'WangAtlas_V3v_left'
                   'WangAtlas_hV4_left'
                   'WangAtlas_VO1_left'
                   'lVOTRC' };

% rmroiCell_VOTRC = rmroiCell(:,6,1:2);
% list_roiNames   = {'lVOTRC' };

list_dtNames     = {'Words_English','Words_Hebrew'};
list_rmNames     = {'retModel-Words_English-css-fFit.mat'
                    'retModel-Words_Hebrew-css-fFit.mat' };

% Launch the function
fname = 'Coverage_EngCB_';  %'Fig1_'; % '' for not saving
fname = '';
ve    = 0.05;
[RF_mean, RF_individuals,empties] = figFunction_coverage_maxProfile_group(cr,list_subInds, ...
                                      'flip',false, ...
                                      'bootcontour', false, ...
                                      'rmroiCell',rmroiCell_noIPS, ...
                                      'list_roiNames', list_roiNames, ...
                                      'list_dtNames', list_dtNames, ...
                                      'list_rmNames', list_rmNames, ...
                                      'sizedegs',7,...
                                      'minvarexp', ve, ...
                                      'numboots',25, ...
                                      'fname', fname, ...
                                      'vers',['v02_' whatFit 'fit'],...
                                      'invisible',true);
%}

%{
% PLOT the mean values per every ROI separately, 3 plots per ROI
% Plot them
mrvNewGraphWin('EngvsHebFOV','wide',true);

ha = tight_subplot(3,6,[.01 .03],[.1 .01],[.01 .01])
for nr=1:length(list_roiNames)
    roiName = list_roiNames{nr};
    ENG     = RF_individuals{nr,1};
    HEB     = RF_individuals{nr,2};
    EMPENG  = empties{nr,1};
    EMPHEB  = empties{nr,2};
    if (isempty(EMPENG) && isempty(EMPHEB))
        alleng = ENG;
        allheb = HEB;
    else
        switch nr
            case 1
                % Some subjects are missing, fix it
                engind = [1,2,3,4,5,6,7,8,9,10,11,12,13];
                hebind = [  1,2,3,  4,5,6,7, 8, 9,10,11];
                engind = [  2,3,4,  6,7,8,9,10,11,12,13];

                % Make paired subjects
                alleng = ENG(:,:,engind);
                allheb = HEB;
                assert(size(ENG,3),size(HEB,3))
            case 2
                % Some subjects are missing, fix it
                engind = [1,2,3,4,5,6,7,8,9,10,11,12,13];
                hebind = [  1,2,3,4,5,6,7,8, 9,10,11,12];
                engind = [  2,3,4,5,6,7,8,9,10,11,12,13];

                % Make paired subjects
                alleng = ENG(:,:,engind);
                allheb = HEB;
                assert(size(ENG,3),size(HEB,3))
            case {3,4,5}
                % Some subjects are missing, fix it
                engind = [1,2,3,4,5,6,7,8,9,10,11,12,13];
                hebind = [  1,2,    3,4,5,6, 7, 8, 9,10];
                engind = [  2,3,    6,7,8,9,10,11,12,13];

                % Make paired subjects
                alleng = ENG(:,:,engind);
                allheb = HEB;
                assert(size(ENG,3),size(HEB,3))
            case 6
                % Some subjects are missing, fix it
                engind = [1,2,3,4,5,6,7,8,9,10,11,12,13];
                hebind = [  1,2,    3,  4,5, 6, 7, 8, 9];
                engind = [  2,3,    6,  7,8, 9,10,11,12];

                % Make paired subjects
                alleng = ENG(:,:,engind);
                allheb = HEB;
                assert(size(ENG,3),size(HEB,3))
            otherwise
                error('This case does not exist')

        end
    end





    % Calculate measures
    mval   = mean(alleng-allheb, 3);
    stdval = std(alleng-allheb, [],3);
    Cd     = zeros(128,128);
    for ii=1:128;for jj=1:128
        Cd(ii,jj)=computeCohen_d(alleng(ii,jj,:),allheb(ii,jj,:),'paired');
    end;end





    axes(ha(nr));
    imagesc(mval);axis equal;colormap(jet);colorbar;grid
    title(sprintf('%s >> Mean of diffs [EngFOV-HebFOV]',strrep(strrep(roiName,'WangAtlas_',''),'_','\_')))
    xlim([1,128]);ylim([1,128])
    xticks([1,64,128]); yticks([1,64,128])
    xticklabels([-7,0,7]); yticklabels([-7,0,7])
    caxis([-.2,.2])
    xlabel('Degs'); ylabel('Degs')

    axes(ha(nr+1*6));
    imagesc(stdval);axis equal;colormap(jet);colorbar;grid
    title(sprintf('%s >> SD of diffs [EngFOV-HebFOV]',strrep(strrep(roiName,'WangAtlas_',''),'_','\_')))
    xlim([1,128]);ylim([1,128])
    xticks([1,64,128]); yticks([1,64,128])
    xticklabels([-7,0,7]); yticklabels([-7,0,7])
    caxis([0,.45])
    xlabel('Degs'); ylabel('Degs')

    axes(ha(nr+2*6));
    imagesc(Cd);axis equal;colormap(jet);colorbar;grid
    title(sprintf("%s >> Cohen's d [EngFOV-HebFOV]",strrep(strrep(roiName,'WangAtlas_',''),'_','\_')))
    xlim([1,128]);ylim([1,128])
    xticks([1,64,128]); yticks([1,64,128])
    xticklabels([-7,0,7]); yticklabels([-7,0,7])
    caxis([-1,1])
    xlabel('Degs'); ylabel('Degs')
end
set(ha(1:12),'XTickLabel',''); set(ha,'YTickLabel','')
set(ha(1:12),'XTick',''); set(ha,'YTick','')
KK=ha(1).XLabel;
set(ha(1:12),'XLabel',KK); set(ha,'YLabel',KK)


KK.String='Degs';
set(ha(13:end),'XTickLabel',[-7,0,7]); set(ha([1,7,12]),'YTickLabel',[-7,0,7])
set(ha(13:end),'XTick',[1,64,128]); set(ha([1,7,12]),'YTick',[1,64,128])

%}


%{
% TESTS JUST VOTRC
rmroiCell_VOTRC = rmroiCell(:,6,1:2);
list_roiNames   ={'lVOTRC' };
list_dtNames     = {'Words_English','Words_Hebrew'};
list_rmNames     = {'retModel-Words_English-css-fFit.mat'
                        'retModel-Words_Hebrew-css-fFit.mat' };
% Launch the function
fname = 'Coverage_EngHeb_';  %'Fig1_'; % '' for not saving
[RF_mean, RF_individuals,empties] = figFunction_coverage_maxProfile_group(cr,list_subInds, 'flip',true, ...
                                      'rmroiCell',rmroiCell_VOTRC, ...
                                      'list_roiNames', list_roiNames, ...
                                      'list_dtNames', list_dtNames, ...
                                      'list_rmNames', list_rmNames, ...
                                      'sizedegs',7,...
                                      'minvarexp', 0.2, ...
                                      'fname', '', ...
                                      'vers',['v01_' whatFit 'fit'],...
                                      'invisible',true);
%}

%{
% PLOT THEM FOR VOTRC, BUT LEAVING ONE SUBJECT AT A TIME
% engind = [  2,3,    6,  7,8, 9,10,11,12];
ALLeng=RF_individuals{1}(:,:,engind);
ALLheb=RF_individuals{2};
% mrvNewGraphWin('EngvsHebFOVreplacement',[],true);
ha = tight_subplot(3,9,[.01 .03],[.1 .01],[.01 .01])
for kk=1:9
    % Remove 1 each time and create same plots with the remaining one
    alleng = ALLeng(:,:,[1:kk-1,kk+1:9]);
    allheb = ALLheb(:,:,[1:kk-1,kk+1:9]);

    % Calculate measures
    mval   = mean(alleng - allheb, 3);
    stdval = std(alleng  - allheb, [],3);
    Cd     = zeros(128,128);
    for ii=1:128;for jj=1:128
            Cd(ii,jj)=computeCohen_d(alleng(ii,jj,:),allheb(ii,jj,:),'paired');
    end;end

    % PLOT
    % subplot(3,9,kk) % (kk*3)-2)
    axes(ha(kk));
    imagesc(mval);axis equal;colormap(jet);colorbar;grid
    title('Mean diffs [Eng-Heb]')
    xlim([1,128]);ylim([1,128])
    % xticks([1,64,128]); yticks([1,64,128])
    % xticklabels([-7,0,7]); yticklabels([-7,0,7])
    % xlabel('Degs'); ylabel('Degs')

    % subplot(3,9,kk+9) % (kk*3)-1)
    axes(ha(kk+9));
    imagesc(stdval);axis equal;colormap(jet);colorbar;grid
    title('SD diffs [Eng-Heb]')
    xlim([1,128]);ylim([1,128])
    % xticks([1,64,128]); yticks([1,64,128])
    % xticklabels([-7,0,7]); yticklabels([-7,0,7])
    % xlabel('Degs'); ylabel('Degs')

    % subplot(3,9,kk+2*9) % kk*3)
    axes(ha(kk+2*9));
    imagesc(Cd);axis equal;colormap(jet);colorbar;grid
    title("Cohen's d")
    xlim([1,128]);ylim([1,128])
    % xticks([1,64,128]); yticks([1,64,128])
    % xticklabels([-7,0,7]); yticklabels([-7,0,7])
    % xlabel('Degs'); ylabel('Degs')
end
set(ha(1:18),'XTickLabel',''); set(ha,'YTickLabel','')
set(ha(1:18),'XTick',''); set(ha,'YTick','')
KK=ha(1).XLabel;
set(ha(1:18),'XLabel',KK); set(ha,'YLabel',KK)


KK.String='Degs';
set(ha(19:end),'XTickLabel',[-7,0,7]); set(ha([1,10,19]),'YTickLabel',[-7,0,7])
set(ha(19:end),'XTick',[1,64,128]); set(ha([1,10,19]),'YTick',[1,64,128])
set(ha(19:end),'XLabel',KK); set(ha([1,10,19]),'YLabel',KK)
%}

% PLOT
%{
mrvNewGraphWin('CrossValEngvsHebFOV',[],true);

subplot(1,3,1)
% subplot(1,2,1)
imagesc(mval);axis equal;colormap(parula);colorbar;grid
title('MEAN: mean of 50 crossvalid. [EngFOV-HebFOV]')
xlim([1,128]);ylim([1,128])
xticks([1,64,128]); yticks([1,64,128])
xticklabels([-7,0,7]); yticklabels([-7,0,7])
xlabel('Degs'); ylabel('Degs')
% subplot(1,2,2)
% surf(mval)




subplot(1,3,2)
imagesc(stdval);axis equal;colormap(parula);colorbar;grid
title('SD: mean of 50 crossvals [EngFOV-HebFOV]')
xlim([1,128]);ylim([1,128])
xticks([1,64,128]); yticks([1,64,128])
xticklabels([-7,0,7]); yticklabels([-7,0,7])
xlabel('Degs'); ylabel('Degs')

% subplot(1,3,3)
%}


%}


rmroiCell_noIPS = rmroiCell(:,1:6,[1,3]);
list_roiNames   = {'WangAtlas_V1v_left'
                   'WangAtlas_V2v_left'
                   'WangAtlas_V3v_left'
                   'WangAtlas_hV4_left'
                   'WangAtlas_VO1_left'
                   'lVOTRC' };

% rmroiCell_VOTRC = rmroiCell(:,6,1:2);
% list_roiNames   = {'lVOTRC' };

list_dtNames     = {'Words_English','Words_Hebrew'};
list_rmNames     = {'retModel-Words_English-css-fFit.mat'
                    'retModel-Words_Hebrew-css-fFit.mat' };

list_dtNames     = {'Words_English','CB'};
list_rmNames     = {'retModel-Words_English-css-fFit.mat'
                    'retModel-Checkers-css-fFit.mat' };

% Launch the function
ves    = [0.05,0.2];
cr.defaults.covfig.vfc.eccthresh = [0.2000 7];
for ve=ves
    fname = ['Coverage_EngCB_13subs-VE' num2str(100*ve)];
    [RF_mean, RF_individuals,empties] = figFunction_coverage_maxProfile_group(...
                                      cr,list_subInds, ...
                                      'flip',false, ...
                                      'bootcontour', false, ...
                                      'rmroiCell',rmroiCell_noIPS, ...
                                      'list_roiNames', list_roiNames, ...
                                      'list_dtNames', list_dtNames, ...
                                      'list_rmNames', list_rmNames, ...
                                      'sizedegs',7,...
                                      'minvarexp', ve, ...
                                      'numboots',25, ...
                                      'fname', fname, ...
                                      'vers',['v02_' whatFit 'fit'],...
                                      'invisible',true);
end
    allsubnames = {'Sub1','Sub2','Sub3','Sub4','Sub5','Sub6','Sub7',...
                           'Sub8','Sub9','Sub10','Sub11','Sub12','Sub13'};
    % Filter results, not all subjects and depending on VE
    RFs = cell(size(RF_individuals));
    switch ve
        case 0.2
            % ROI 1
            % Some subjects are missing, fix it
             hebind{1}   = [1:13];
             engind{1}   = [1:13];
             subnames{1} = allsubnames;
             hebind{2}=hebind{1};engind{2}=engind{1};subnames{2}=subnames{1};
             hebind{3}=hebind{1};engind{3}=engind{1};subnames{3}=subnames{1};
             hebind{4}=hebind{1};engind{4}=engind{1};subnames{4}=subnames{1};
             hebind{5}=hebind{1};engind{5}=engind{1};subnames{5}=subnames{1};

%             hebind{1}   = [  1,2,3,  4,5,6,7, 8, 9,10,11];
%             engind{1}   = [  2,3,4,  6,7,8,9,10,11,12,13];
%             subnames{1} = allsubnames(engind{1});
%
%
%             % ROI 2
%             % Some subjects are missing, fix it
%             hebind{2}   = [  1,2,3,4,5,6,7,8, 9,10,11,12];
%             engind{2}   = [  2,3,4,5,6,7,8,9,10,11,12,13];
%             subnames{2} = allsubnames(engind{2});
%
%             % ROI {3,4,5}
%             % Some subjects are missing, fix it
%             hebind{3}   = [  1,2,    3,4,5,6, 7, 8, 9,10];
%             engind{3}   = [  2,3,    6,7,8,9,10,11,12,13];
%             subnames{3} = allsubnames(engind{3});
%             hebind{4}=hebind{3};engind{4}=engind{3};subnames{4}=subnames{3};
%             hebind{5}=hebind{3};engind{5}=engind{3};subnames{5}=subnames{3};
%
%             % ROI 6
%             hebind{6}   = [  1,2,    3,  4,5, 6, 7, 8, 9];
%             engind{6}   = [  2,3,    6,  7,8, 9,10,11,12];
%             subnames{6} = {'Sub2','Sub3','Sub6','Sub8','Sub9','Sub10',...
%                 'Sub11','Sub12','Sub13'};

            hebind{6}   = [  1,2,3,4,5,6,8,9,10,11,12,13];
            engind{6}   = [  1,2,3,4,5,6,7,8,9,10,11,12];
            subnames{6} = {'Sub1','Sub2','Sub3','Sub4','Sub5','Sub6',...
                           'Sub8','Sub9','Sub10','Sub11','Sub12','Sub13'};

        case 0.05
            % ROIs 1 to 5
            hebind{1}   = [1:13];
            engind{1}   = [1:13];
            subnames{1} = {'Sub1','Sub2','Sub3','Sub4','Sub5','Sub6','Sub7',...
                'Sub8','Sub9','Sub10','Sub11','Sub12','Sub13'};
            hebind{2}=hebind{1};engind{2}=engind{1};subnames{2}=subnames{1};
            hebind{3}=hebind{1};engind{3}=engind{1};subnames{3}=subnames{1};
            hebind{4}=hebind{1};engind{4}=engind{1};subnames{4}=subnames{1};
            hebind{5}=hebind{1};engind{5}=engind{1};subnames{5}=subnames{1};
            hebind{6}=hebind{1};engind{6}=engind{1};subnames{6}=subnames{1};
            % ROI 6
%             hebind{6}   = [1:12];
%             engind{6}   = [1:3,5:13];
%             subnames{6} = {'Sub1','Sub2','Sub3','Sub5','Sub6','Sub7','Sub8',...
%                            'Sub9','Sub10','Sub11','Sub12','Sub13'};
    end
    for ii=1:6
        % Eng
        tmp       = RF_individuals(ii,1);
        RFs(ii,1) = {tmp{1}(:,:,engind{ii})};
        % Heb
        tmp       = RF_individuals(ii,2);
        RFs(ii,2) = {tmp{1}(:,:,hebind{ii})};
    end


    % BOOTSTRAPPING
    bootstrapping = true;
    for numr=1:6
        MVALS = zeros(128,128,50);
        SVALS = zeros(128,128,50);
        DVALS = zeros(128,128,50);

        ALLeng = RFs{numr,1};
        ALLheb = RFs{numr,2};

        if bootstrapping
            for kk=1:50
                % Remove 1 each time and create same plots with the remaining one
                % randReplacement = datasample(1:9,9);
                % Control the ranzomization with rng so that this is reproducible
                % In checkers, only survive 8
                rng(kk)
                randReplacement = datasample(1:size(ALLeng,3),size(ALLeng,3));

                alleng          = ALLeng(:,:,randReplacement);
                allheb          = ALLheb(:,:,randReplacement);

                % Calculate measures
                tmval   = mean(alleng - allheb, 3);
                tstdval = std(alleng  - allheb, [],3);
                tCd     = zeros(128,128);
                for ii=1:128; for jj=1:128
                    tCd(ii,jj)=computeCohen_d(alleng(ii,jj,:),allheb(ii,jj,:),'paired');
                end; end

                % Accummlate it
                MVALS(:,:,kk) = tmval;
                SVALS(:,:,kk) = tstdval;
                DVALS(:,:,kk) = tCd;
            end
            % obtain means again
            mval{numr}   = mean(MVALS,3);
            stdval{numr} = mean(SVALS,3);
            Cd{numr}     = mean(DVALS,3);
        else
            % Calculate measures
            mval{numr}   = mean(ALLeng - ALLheb, 3);
            stdval{numr} = std(ALLeng  - ALLheb, [],3);
            Cd{numr}     = zeros(128,128);
            for ii=1:128; for jj=1:128
                    Cd{numr}(ii,jj)=computeCohen_d(ALLeng(ii,jj,:),ALLheb(ii,jj,:),'paired');
            end;end
        end
    end  % numr

    % PLOTS
% end  % End VE

% for ve=ves
    % Plot d' in all ROIs
    mrvNewGraphWin('alldprimes','wide',true);
    ha = tight_subplot(1,6,[.01 .03],[.1 .01],[.01 .01]);
    for nr=1:length(list_roiNames)
        roiName = list_roiNames{nr};

        %{
        axes(ha(nr));
        imagesc(mval);axis equal;colormap(jet);colorbar;grid
        title(sprintf('%s >> Mean of diffs [EngFOV-HebFOV]',strrep(strrep(roiName,'WangAtlas_',''),'_','\_')))
        xlim([1,128]);ylim([1,128])
        xticks([1,64,128]); yticks([1,64,128])
        xticklabels([-7,0,7]); yticklabels([-7,0,7])
        caxis([-.2,.2])
        xlabel('Degs'); ylabel('Degs')

        axes(ha(nr+1*6));
        imagesc(stdval);axis equal;colormap(jet);colorbar;grid
        title(sprintf('%s >> SD of diffs [EngFOV-HebFOV]',strrep(strrep(roiName,'WangAtlas_',''),'_','\_')))
        xlim([1,128]);ylim([1,128])
        xticks([1,64,128]); yticks([1,64,128])
        xticklabels([-7,0,7]); yticklabels([-7,0,7])
        caxis([0,.45])
        xlabel('Degs'); ylabel('Degs')
        %}

        axes(ha(nr));
        imagesc(Cd{nr});axis equal;colormap(jet);colorbar;grid
        title(sprintf("%s >> Cohen's d [EngFOV-HebFOV]",strrep(strrep(roiName,'WangAtlas_',''),'_','\_')))
        xlim([1,128]);ylim([1,128])
        xticks([1,64,128]); yticks([1,64,128])
        xticklabels([-7,0,7]); yticklabels([-7,0,7])
        % caxis([round(min(Cd{nr}(:)),2),round(max(Cd{nr}(:)),2)])
        caxis([-1.5,1.5])
        xlabel('Degs'); ylabel('Degs')
    end
    set(ha(2:6),'YTickLabel',''); set(ha(2:6),'YTick','');
    set(ha(1:end),'XTickLabel',[-7,0,7]); set(ha([1]),'YTickLabel',[-7,0,7])
    set(ha(1:end),'XTick',[1,64,128]); set(ha([1]),'YTick',[1,64,128])

    titlefile  = ['alldprimes_Eng-Heb-' num2str(size(ALLeng,3)) ...
                  'subs-VE' num2str(100*ve)];
    saveas(gcf, fullfile(cr.dirs.FIGPNG,[titlefile '.png']), 'png')
    saveas(gcf, fullfile(cr.dirs.FIGSVG,[titlefile '.svg']), 'svg')

    % saveas(gcf, fullfile(sdRP,'DATA','figures','png',[titlefile '.png']), 'png')
    % saveas(gcf, fullfile(sdRP,'DATA','figures','svg',[titlefile '.svg']), 'svg')
    close all









    % Plot the mesh and the d' in the VOTRC
    for nnrr=6
        mrvNewGraphWin('CrossValEngvsHebFOV',[],true);

        cdlims = [-1.5,1.5];

        subplot(1,2,1)
        imagesc(Cd{6});axis equal;colormap(parula);colorbar;grid
        title("Cohen's d: mean of 50 crossvals [Eng-Heb]")
        xlim([1,128]);ylim([1,128])
        xticks([1,64,128]); yticks([1,64,128])
        % caxis([-0.7,1.7])
        caxis(cdlims)
        xticklabels([-7,0,7]); yticklabels([-7,0,7])
        xlabel('Degs'); ylabel('Degs')


        subplot(1,2,2)

        [X,Y] = meshgrid(1:128,1:128);
        XX = ((X-64)/64)*7;
        YY = ((Y-64)/64)*7;
        YY = flipud(YY);
        surf(XX,YY,Cd{6});
        xlabel('X (degs)'); ylabel('Y (degs)')
        zlabel("Cohen's d")
        xlim([-7,7])
        ylim([-7,7])
        zlim(cdlims)
        xticks([-7,-5,-3,-1,1,3,5,7])
        yticks([-7,-5,-3,-1,1,3,5,7])
        xticklabels({'-7','-5','-3','-1','1','3','5','7'})
        yticklabels({''})
        set(gca,'FontSize',18)


        titlefile  = ['Meshdprime_' list_roiNames{nnrr} ...
                      '_Eng-Heb-' num2str(size(ALLeng,3)) ...
                      'subs-VE' num2str(100*ve)];
        saveas(gcf, fullfile(cr.dirs.FIGPNG,[titlefile '.png']), 'png')
        saveas(gcf, fullfile(cr.dirs.FIGSVG,[titlefile '.svg']), 'svg')
        % saveas(gcf, fullfile(sdRP,'DATA','figures','png',[titlefile '.png']), 'png')
        % saveas(gcf, fullfile(sdRP,'DATA','figures','fig',[titlefile '.fig']), 'fig')
        % saveas(gcf, fullfile(sdRP,'DATA','figures','svg',[titlefile '.svg']), 'svg')
        close all
    end

    % Plot individual subject differences in VOTRC
    mrvNewGraphWin('CrossValEngvsHebFOV','wide',true);
    switch ve
        case 0.2
            % co=20%
            ha = tight_subplot(1,size(ALLeng,3),[.01 .03],[.1 .01],[.01 .01]);
        case 0.05
            % co=5%
            ha = tight_subplot(2,size(ALLeng,3)/2,[.01 .03],[.1 .01],[.01 .01]);
    end

    ALLeng = RFs{6,1};
    ALLheb = RFs{6,2};

    for nn=1:size(ALLeng,3)
        ieng = ALLeng(:,:,nn);
        iheb = ALLheb(:,:,nn);
        axes(ha(nn));
        imagesc(ieng-iheb);axis equal;colormap(jet);colorbar;grid
        caxis([-1,1])
        % title(sprintf('Sub ind %i',nn))
        title(subnames{6}{nn})
        xlim([1,128]);ylim([1,128])
        xticks([1,64,128]); yticks([1,64,128])
        xticklabels([-7,0,7]); yticklabels([-7,0,7])

    end

    ha    = xlabel('Degs');
    ha(1) = ylabel('Degs');
    titlefile  = ['IndividualSubjectEng-Heb-' num2str(size(ALLeng,3)) ...
                                                     'subs-VE' num2str(100*ve)];
    saveas(gcf, fullfile(cr.dirs.FIGPNG,[titlefile '.png']), 'png')
    saveas(gcf, fullfile(cr.dirs.FIGSVG,[titlefile '.svg']), 'svg')
    % saveas(gcf, fullfile(sdRP,'DATA','figures','png',[titlefile '.png']), 'png')
    close all
% end

%% FIGURE S4: WE_CB
% Order is WE_CB
list_roiNames16 = list_roiNames(1:8);
list_rmDescripts = strrep({'Words_English','Checkers'},'_','\_');

% Obtain equally thresholded voxels to scatterplot
ves=[0.20,0.05];
for ve=ves
    [R,C_data,cr]=crThreshGetSameVoxel(cr,...
                                   rmroiCell,...
                                   list_subInds,...
                                   list_roiNames16,...
                                   'cothres', ve,...
                                   'fieldrange', 7);

%{
    fprintf('------------------------------------\n')
    fprintf('         %s - %s (R2:%g)        \n', list_rmDescripts{2},list_rmDescripts{1},cothresh)
    fprintf('------------------------------------\n\n')
    for ii=1:8
        fprintf('%s(N=%i)\n',strrep(list_roiNames16{ii},'WangAtlas_',''),length(R.X_rm2{ii}))
        [H P CI] = ttest(R.X_rm2{ii}-R.X_rm1{ii});
        fprintf('(X) P:%g, CI: [%g %g]\n',P,CI(1),CI(2))

        [H P CI] = ttest(R.Y_rm2{ii}-R.Y_rm1{ii});
        fprintf('(Y) P:%g, CI: [%g %g]\n\n',P,CI(1),CI(2))
    end
%}

    fname = ['scatterplot_eccentricity_WordEngVsCB_VE' num2str(100*ve) '_8ROIs_13subs_' whatFit 'Fit_v02'];
    % fname = '';
    crCreateScatterplot(R,C_data,cr,...
                        list_subInds,...
                        list_roiNames16,...
                        list_rmDescripts,...
                        'ecc', ...  % 'co'
                        fontsize, ...
                        fname);
    fname = ['scatterplot_varianceExplained_WordEngVsCB_VE' num2str(100*ve) '_8ROIs_13subs_' whatFit 'Fit_v02'];
    crCreateScatterplot(R,C_data,cr,...
                        list_subInds,...
                        list_roiNames,...
                        list_rmDescripts,...
                        'co', ...
                        fontsize,...
                        fname);
end
%{
crCreateScatterplot(R,C_data,cr,...
                    list_subInds,...
                    list_roiNames16,...
                    list_rmDescripts,...
                    'x0', ...  % 'co'
                    fontsize, ...
                    fname);

crCreateScatterplot(R,C_data,cr,...
                    list_subInds,...
                    list_roiNames16,...
                    list_rmDescripts,...
                    'y0', ...  % 'co'
                    fontsize, ...
                    fname);

 %}


%% FIGURE 6: (B) Variance Explained: Scatterplot: wordLarge-WordSmall
% Uses the same voxel calculations from the previous plot
% If doubt or this is moved, calculate it again here

% Plot it
fname = ['scatterplot_varianceExplained_WordEngVsWordHeb_6ROIs_13subs_' whatFit 'Fit_v02'];
fname = '';
crCreateScatterplot(R,C_data,cr,...
                    list_subInds,...
                    list_roiNames16,...
                    list_rmDescripts,...
                    'co', ...  % 'co', 'ecc'
                    fname);

% Group COVERAGE plots, take all subjects from list_subInds
rmroiCell_VOTRC  = rmroiCell(:,6,1:2);
list_roiNames    = {'lVOTRC' };
list_dtNames     = {'Words_English','Checkers'};
list_rmNames     = {'retModel-Words_English-css-fFit.mat'
                    'retModel-Checkers-css-fFit.mat' };

% Launch the function
fname = 'Coverage_EngCB_new_';  %'Fig1_'; % '' for not saving
% fname = '';
[RF_mean, RF_individuals,empties] = figFunction_coverage_maxProfile_group(cr,list_subInds, ...
                                      'flip',false, ...
                                      'bootcontour', false, ...
                                      'rmroiCell',rmroiCell_VOTRC, ...
                                      'list_roiNames', list_roiNames, ...
                                      'list_dtNames', list_dtNames, ...
                                      'list_rmNames', list_rmNames, ...
                                      'sizedegs',7,...
                                      'minvarexp', 0.05, ...
                                      'numboots',25, ...
                                      'fname', fname, ...
                                      'vers',['v02_' whatFit 'fit'],...
                                      'invisible',false);

% PLOT the mean values per every ROI separately, 3 plots per ROI
%{
% Plot them
mrvNewGraphWin('EngvsHebFOV','wide',true);

ha = tight_subplot(3,6,[.01 .03],[.1 .01],[.01 .01])
for nr=1:length(list_roiNames)
    roiName = list_roiNames{nr};
    ENG     = RF_individuals{nr,1};
    HEB     = RF_individuals{nr,2};
    EMPENG  = empties{nr,1};
    EMPHEB  = empties{nr,2};
    if (isempty(EMPENG) && isempty(EMPHEB))
        alleng = ENG;
        allheb = HEB;
    else
        switch nr
            case 1
                % Some subjects are missing, fix it
                engind = [1,2,3,4,5,6,7,8,9,10,11,12,13];
                hebind = [  1,2,3,  4,5,6,7, 8, 9,10,11];
                engind = [  2,3,4,  6,7,8,9,10,11,12,13];

                % Make paired subjects
                alleng = ENG(:,:,engind);
                allheb = HEB;
                assert(size(ENG,3),size(HEB,3))
            case 2
                % Some subjects are missing, fix it
                engind = [1,2,3,4,5,6,7,8,9,10,11,12,13];
                hebind = [  1,2,3,4,5,6,7,8, 9,10,11,12];
                engind = [  2,3,4,5,6,7,8,9,10,11,12,13];

                % Make paired subjects
                alleng = ENG(:,:,engind);
                allheb = HEB;
                assert(size(ENG,3),size(HEB,3))
            case {3,4,5}
                % Some subjects are missing, fix it
                engind = [1,2,3,4,5,6,7,8,9,10,11,12,13];
                hebind = [  1,2,    3,4,5,6, 7, 8, 9,10];
                engind = [  2,3,    6,7,8,9,10,11,12,13];

                % Make paired subjects
                alleng = ENG(:,:,engind);
                allheb = HEB;
                assert(size(ENG,3),size(HEB,3))
            case 6
                % Some subjects are missing, fix it
                engind = [1,2,3,4,5,6,7,8,9,10,11,12,13];
                hebind = [  1,2,    3,  4,5, 6, 7, 8, 9];
                engind = [  2,3,    6,  7,8, 9,10,11,12];

                % Make paired subjects
                alleng = ENG(:,:,engind);
                allheb = HEB;
                assert(size(ENG,3),size(HEB,3))
            otherwise
                error('This case does not exist')

        end
    end

    % Calculate measures
    mval   = mean(alleng-allheb, 3);
    stdval = std(alleng-allheb, [],3);
    Cd     = zeros(128,128);
    for ii=1:128;for jj=1:128
        Cd(ii,jj)=computeCohen_d(alleng(ii,jj,:),allheb(ii,jj,:),'paired');
    end;end





    axes(ha(nr));
    imagesc(mval);axis equal;colormap(jet);colorbar;grid
    title(sprintf('%s >> Mean of diffs [EngFOV-HebFOV]',strrep(strrep(roiName,'WangAtlas_',''),'_','\_')))
    xlim([1,128]);ylim([1,128])
    xticks([1,64,128]); yticks([1,64,128])
    xticklabels([-7,0,7]); yticklabels([-7,0,7])
    caxis([-.2,.2])
    xlabel('Degs'); ylabel('Degs')

    axes(ha(nr+1*6));
    imagesc(stdval);axis equal;colormap(jet);colorbar;grid
    title(sprintf('%s >> SD of diffs [EngFOV-HebFOV]',strrep(strrep(roiName,'WangAtlas_',''),'_','\_')))
    xlim([1,128]);ylim([1,128])
    xticks([1,64,128]); yticks([1,64,128])
    xticklabels([-7,0,7]); yticklabels([-7,0,7])
    caxis([0,.45])
    xlabel('Degs'); ylabel('Degs')

    axes(ha(nr+2*6));
    imagesc(Cd);axis equal;colormap(jet);colorbar;grid
    title(sprintf("%s >> Cohen's d [EngFOV-HebFOV]",strrep(strrep(roiName,'WangAtlas_',''),'_','\_')))
    xlim([1,128]);ylim([1,128])
    xticks([1,64,128]); yticks([1,64,128])
    xticklabels([-7,0,7]); yticklabels([-7,0,7])
    caxis([-1,1])
    xlabel('Degs'); ylabel('Degs')
end
set(ha(1:12),'XTickLabel',''); set(ha,'YTickLabel','')
set(ha(1:12),'XTick',''); set(ha,'YTick','')
KK=ha(1).XLabel;
set(ha(1:12),'XLabel',KK); set(ha,'YLabel',KK)


KK.String='Degs';
set(ha(13:end),'XTickLabel',[-7,0,7]); set(ha([1,7,12]),'YTickLabel',[-7,0,7])
set(ha(13:end),'XTick',[1,64,128]); set(ha([1,7,12]),'YTick',[1,64,128])

%}

% TESTS JUST VOTRC
%{
rmroiCell_VOTRC = rmroiCell(:,6,1:2);
list_roiNames   ={'lVOTRC' };
list_dtNames     = {'Words_English','Words_Hebrew'};
list_rmNames     = {'retModel-Words_English-css-fFit.mat'
                        'retModel-Words_Hebrew-css-fFit.mat' };
% Launch the function
fname = 'Coverage_EngHeb_';  %'Fig1_'; % '' for not saving
[RF_mean, RF_individuals,empties] = figFunction_coverage_maxProfile_group(cr,list_subInds, 'flip',true, ...
                                      'rmroiCell',rmroiCell_VOTRC, ...
                                      'list_roiNames', list_roiNames, ...
                                      'list_dtNames', list_dtNames, ...
                                      'list_rmNames', list_rmNames, ...
                                      'sizedegs',7,...
                                      'minvarexp', 0.2, ...
                                      'fname', '', ...
                                      'vers',['v01_' whatFit 'fit'],...
                                      'invisible',true);
%}

% PLOT THEM FOR VOTRC, BUT LEAVING ONE SUBJECT AT A TIME
%{
% engind = [  2,3,    6,  7,8, 9,10,11,12];
ALLeng=RF_individuals{1}(:,:,engind);
ALLheb=RF_individuals{2};
% mrvNewGraphWin('EngvsHebFOVreplacement',[],true);
ha = tight_subplot(3,9,[.01 .03],[.1 .01],[.01 .01])
for kk=1:9
    % Remove 1 each time and create same plots with the remaining one
    alleng = ALLeng(:,:,[1:kk-1,kk+1:9]);
    allheb = ALLheb(:,:,[1:kk-1,kk+1:9]);

    % Calculate measures
    mval   = mean(alleng - allheb, 3);
    stdval = std(alleng  - allheb, [],3);
    Cd     = zeros(128,128);
    for ii=1:128;for jj=1:128
            Cd(ii,jj)=computeCohen_d(alleng(ii,jj,:),allheb(ii,jj,:),'paired');
    end;end

    % PLOT
    % subplot(3,9,kk) % (kk*3)-2)
    axes(ha(kk));
    imagesc(mval);axis equal;colormap(jet);colorbar;grid
    title('Mean diffs [Eng-Heb]')
    xlim([1,128]);ylim([1,128])
    % xticks([1,64,128]); yticks([1,64,128])
    % xticklabels([-7,0,7]); yticklabels([-7,0,7])
    % xlabel('Degs'); ylabel('Degs')

    % subplot(3,9,kk+9) % (kk*3)-1)
    axes(ha(kk+9));
    imagesc(stdval);axis equal;colormap(jet);colorbar;grid
    title('SD diffs [Eng-Heb]')
    xlim([1,128]);ylim([1,128])
    % xticks([1,64,128]); yticks([1,64,128])
    % xticklabels([-7,0,7]); yticklabels([-7,0,7])
    % xlabel('Degs'); ylabel('Degs')

    % subplot(3,9,kk+2*9) % kk*3)
    axes(ha(kk+2*9));
    imagesc(Cd);axis equal;colormap(jet);colorbar;grid
    title("Cohen's d")
    xlim([1,128]);ylim([1,128])
    % xticks([1,64,128]); yticks([1,64,128])
    % xticklabels([-7,0,7]); yticklabels([-7,0,7])
    % xlabel('Degs'); ylabel('Degs')
end
set(ha(1:18),'XTickLabel',''); set(ha,'YTickLabel','')
set(ha(1:18),'XTick',''); set(ha,'YTick','')
KK=ha(1).XLabel;
set(ha(1:18),'XLabel',KK); set(ha,'YLabel',KK)


KK.String='Degs';
set(ha(19:end),'XTickLabel',[-7,0,7]); set(ha([1,10,19]),'YTickLabel',[-7,0,7])
set(ha(19:end),'XTick',[1,64,128]); set(ha([1,10,19]),'YTick',[1,64,128])
set(ha(19:end),'XLabel',KK); set(ha([1,10,19]),'YLabel',KK)

%}



% PLOT THEM FOR VOTRC, DO BOOTSTRAPPING AND AVERAGE IT
% engind = [  2,3,    6,  7,8, 9,10,11,12];
% When using checkers
% old fit
%{
allind = [1,2,3,4,5,6,7,8,9,10,11,12,13];
chcind = [    1,  2,    3,4, 5, 6, 7, 8];
engind = [1,2,3,4,5,6,  7,8, 9,10,11,12];
engind = [    3  ,5,    7,8, 9,10,11,12];
%}
% new fit
allind = [1,2,3,4,5,6,7,8,9,10,11,12,13];
% R2 > 20%
chcind = [    1,  2,    3,4, 5, 6, 7   ];
engind = [1,2,3,4,5,6,  7,8, 9,10,11,12];
engind = [    3  ,5,    7,8, 9,10,12];
% R2 > 10%
chcind = [  1,2,3,4,    5,6, 7, 8, 9,10];
engind = [1,2,3,4,5,6,7,8,9,10,11,12,13];
engind = [  2,3,4,5,    8,9,10,11,12,13];


ALLeng=RF_individuals{1}(:,:,engind);
ALLheb=RF_individuals{2};

MVALS = zeros(128,128,50);
SVALS = zeros(128,128,50);
DVALS = zeros(128,128,50);

for kk=1:100
    % Remove 1 each time and create same plots with the remaining one
    assert(isequal(size(ALLeng,3),size(ALLheb,3)))
    randReplacement = datasample(1:size(ALLheb,3),size(ALLheb,3));

    alleng          = ALLeng(:,:,randReplacement);
    allheb          = ALLheb(:,:,randReplacement);

    % Calculate measures
    mval   = mean(alleng - allheb, 3);
    stdval = std(alleng  - allheb, [],3);
    Cd     = zeros(128,128);
    for ii=1:128;for jj=1:128
            Cd(ii,jj)=computeCohen_d(alleng(ii,jj,:),allheb(ii,jj,:),'paired');
    end;end

    % Accummlate it
    MVALS(:,:,kk) = mval;
    SVALS(:,:,kk) = stdval;
    DVALS(:,:,kk) = Cd;

end

% obtain means again
mval   = mean(MVALS,3);
stdval = mean(SVALS,3);
Cd     = mean(DVALS,3);

% min(Cd(:))
% max(Cd(:))
% min(Cd(:))
% mean(Cd(:))
%
% 100*sum(Cd(:)>0)/length(Cd(:))
% 100*sum(Cd(:)<=0)/length(Cd(:))


mrvNewGraphWin('CrossValEngvsHebFOV','wide',true);

caxismin = min(Cd(:));
caxismax = max(Cd(:));

subplot(1,2,1)
imagesc(Cd);axis equal;colormap(parula);colorbar;grid
title("Cohen's d: mean of 50 crossvals [Eng-CB]")
xlim([1,128]);ylim([1,128])
xticks([1,64,128]); yticks([1,64,128])
caxis([caxismin,caxismax])
xticklabels([-7,0,7]); yticklabels([-7,0,7])
xlabel('Degs'); ylabel('Degs')
% cmap = colormap(parula);
% cmap1 = cmap(1:5:128,:);
% cmap2 = cmap(129:end,:);
% vals  = linspace(0,1,50)';
% lght  = linspace(0,.3,50)';
% cmap1 = [vals,lght,lght];
% cmap2 = [lght,lght,vals];
% cmap=[flip(cmap2);cmap1];
% cmap=[cmap1;cmap2];
% cmap = cmap(129:end,:);
% colormap(cmap);
colorbar





subplot(1,2,2)
[X,Y] = meshgrid(1:128,1:128);
XX = ((X-64)/64)*7;
YY = ((Y-64)/64)*7;
YY = flipud(YY);
surf(XX,YY,Cd);
xlabel('X (degs)'); ylabel('Y (degs)')
zlabel("Cohen's d")


xlim([-7,7])
ylim([-7,7])
zlim([caxismin,caxismax])
xticks([-7,-5,-3,-1,1,3,5,7])
yticks([-7,-5,-3,-1,1,3,5,7])
xticklabels({'-7','-5','-3','-1','1','3','5','7'})
yticklabels({''})
set(gca,'FontSize',18)
caxis([caxismin,caxismax])
% colormap(cmap);
colormap(parula);



%% PREPARE DATA: HEBREW AND CB FROM ISRAEL
% Generate the rmroicell that we will use in all plots in this script
% This will read the results obtained by Rosemary or the re-run in 2021
fontsize = 12;
readExisting = true;
% Do the same with the small and large words
list_subInds     = [31:36 38:44];
list_roiNames = {'WangAtlas_V1v_left'
                 'WangAtlas_V2v_left'
                 'WangAtlas_V3v_left'
                 'WangAtlas_hV4_left'
                 'WangAtlas_VO1_left'
                 'lVOTRC'
                 'WangAtlas_IPS0'
                 'WangAtlas_IPS1'};


whatFit = 'new';  % 'new' | 'Rosemary'
list_dtNames     = {'Words_English','Words_Hebrew','Checkers'};
list_rmDescripts = {'Words_English','Words_Hebrew','Checkers'};
if strcmp(whatFit,'Rosemary')
    list_rmNames     = {'retModel-Words_English-css.mat',...
                        'retModel-Words_Hebrew-css.mat', ...
                        'retModel-Checkers-css.mat'};
else
     list_rmNames     = {'retModel-Words_English-css-fFit.mat',...
                         'retModel-Words_Hebrew-css-fFit.mat', ...
                         'retModel-Checkers-css-fFit.mat'};
end
matname = ['rmroicell_subInds-31to36-38to44_dtNames-ALL_fits-' whatFit '_2022.mat'];


if readExisting
    load(fullfile(sdRP,'DATA',matname),'rmroiCell');

    % If reading checkers
    % rmroiCell = rmroiCell(:,1:6,:);
    % A = load(fullfile(sdRP,'DATA',matname2),'rmroiCell');
    % Overwrite ENGLISH to have HEBREW and CB
    % rmroiCell(:,1:6,1) = A.rmroiCell(:,1:6);
    % Change order so that it is x:Heb, y:CB
    rmroiCell = rmroiCell(:,:,[2,1]);
else
    rmroiCell = ff_rmroiCell(cr,...
                             list_subInds,...
                             list_roiNames,...
                             list_dtNames, ...
                             list_rmNames,...
                             'list_path',cr.bk.list_sessionRet,...
                             'latest_fFit',true, ...
                             'checkYear','2022');
    % Save rmroicell
    sss
end

% Read the generic params for coverage for all subjects
cr.defaults.covfig.vfc = ff_vfcDefault();
cr.defaults.covfig.vfc.list_roiNames    = list_roiNames;
% data types we want to look at
cr.defaults.covfig.vfc.list_dtNames     = list_dtNames;
% names of the rm in each dt
cr.defaults.covfig.vfc.list_rmNames     = list_rmNames;
cr.defaults.covfig.vfc.list_rmDescripts = list_rmDescripts;
% subinds = [31:36 38:44]; % Hebrew
% cr.defaults.covfig.vfc = ff_vfcDefault_Hebrew();


% Order is WE_CB
list_roiNames16 = list_roiNames(1:6);
list_rmDescripts = {'Words_Hebrew','Checkers'};

% Obtain equally thresholded voxels to scatterplot
cothresh = 0.05;
[R,C_data,cr]=crThreshGetSameVoxel(cr,...
                                   rmroiCell,...
                                   list_subInds,...
                                   list_roiNames16,...
                                   'cothres', cothresh,...
                                   'fieldrange', 7);


fprintf('------------------------------------\n')
fprintf('         %s - %s (R2:%g)        \n', list_rmDescripts{2},list_rmDescripts{1},cothresh)
fprintf('------------------------------------\n\n')
for ii=1:6
    fprintf('%s(N=%i)\n',strrep(list_roiNames16{ii},'WangAtlas_',''),length(R.X_rm2{ii}))
    [H P CI] = ttest(R.X_rm2{ii}-R.X_rm1{ii});
    fprintf('(X) P:%g, CI: [%g %g]\n',P,CI(1),CI(2))

    [H P CI] = ttest(R.Y_rm2{ii}-R.Y_rm1{ii});
    fprintf('(Y) P:%g, CI: [%g %g]\n\n',P,CI(1),CI(2))
end


fname = ['scatterplot_eccentricity_WordEngVsWordHeb_6ROIs_13subs_' whatFit 'Fit_v02'];
fname = '';
crCreateScatterplot(R,C_data,cr,...
                    list_subInds,...
                    list_roiNames16,...
                    list_rmDescripts,...
                    'ecc', ...  % 'co'
                    fontsize, ...
                    fname);

crCreateScatterplot(R,C_data,cr,...
                    list_subInds,...
                    list_roiNames16,...
                    list_rmDescripts,...
                    'x0', ...  % 'co'
                    fontsize, ...
                    fname);

crCreateScatterplot(R,C_data,cr,...
                    list_subInds,...
                    list_roiNames16,...
                    list_rmDescripts,...
                    'y0', ...  % 'co'
                    fontsize, ...
                    fname);




% FIGURE 6: (B) Variance Explained: Scatterplot: wordLarge-WordSmall
% Uses the same voxel calculations from the previous plot
% If doubt or this is moved, calculate it again here

% Plot it
fname = ['scatterplot_varianceExplained_WordEngVsWordHeb_6ROIs_13subs_' whatFit 'Fit_v02'];
fname = '';
crCreateScatterplot(R,C_data,cr,...
                    list_subInds,...
                    list_roiNames16,...
                    list_rmDescripts,...
                    'co', ...  % 'co', 'ecc'
                    fname);

% Group COVERAGE plots, take all subjects from list_subInds
rmroiCell_VOTRC  = rmroiCell(:,6,1:2);
list_roiNames    = {'lVOTRC' };
list_dtNames     = {'Words_English','Checkers'};
list_rmNames     = {'retModel-Words_English-css-fFit.mat'
                    'retModel-Checkers-css-fFit.mat' };

% Launch the function
fname = 'Coverage_EngCB_new_';  %'Fig1_'; % '' for not saving
% fname = '';
[RF_mean, RF_individuals,empties] = figFunction_coverage_maxProfile_group(cr,list_subInds, ...
                                      'flip',false, ...
                                      'bootcontour', false, ...
                                      'rmroiCell',rmroiCell_VOTRC, ...
                                      'list_roiNames', list_roiNames, ...
                                      'list_dtNames', list_dtNames, ...
                                      'list_rmNames', list_rmNames, ...
                                      'sizedegs',7,...
                                      'minvarexp', 0.05, ...
                                      'numboots',25, ...
                                      'fname', fname, ...
                                      'vers',['v02_' whatFit 'fit'],...
                                      'invisible',false);

% PLOT the mean values per every ROI separately, 3 plots per ROI
%{
% Plot them
mrvNewGraphWin('EngvsHebFOV','wide',true);

ha = tight_subplot(3,6,[.01 .03],[.1 .01],[.01 .01])
for nr=1:length(list_roiNames)
    roiName = list_roiNames{nr};
    ENG     = RF_individuals{nr,1};
    HEB     = RF_individuals{nr,2};
    EMPENG  = empties{nr,1};
    EMPHEB  = empties{nr,2};
    if (isempty(EMPENG) && isempty(EMPHEB))
        alleng = ENG;
        allheb = HEB;
    else
        switch nr
            case 1
                % Some subjects are missing, fix it
                engind = [1,2,3,4,5,6,7,8,9,10,11,12,13];
                hebind = [  1,2,3,  4,5,6,7, 8, 9,10,11];
                engind = [  2,3,4,  6,7,8,9,10,11,12,13];

                % Make paired subjects
                alleng = ENG(:,:,engind);
                allheb = HEB;
                assert(size(ENG,3),size(HEB,3))
            case 2
                % Some subjects are missing, fix it
                engind = [1,2,3,4,5,6,7,8,9,10,11,12,13];
                hebind = [  1,2,3,4,5,6,7,8, 9,10,11,12];
                engind = [  2,3,4,5,6,7,8,9,10,11,12,13];

                % Make paired subjects
                alleng = ENG(:,:,engind);
                allheb = HEB;
                assert(size(ENG,3),size(HEB,3))
            case {3,4,5}
                % Some subjects are missing, fix it
                engind = [1,2,3,4,5,6,7,8,9,10,11,12,13];
                hebind = [  1,2,    3,4,5,6, 7, 8, 9,10];
                engind = [  2,3,    6,7,8,9,10,11,12,13];

                % Make paired subjects
                alleng = ENG(:,:,engind);
                allheb = HEB;
                assert(size(ENG,3),size(HEB,3))
            case 6
                % Some subjects are missing, fix it
                engind = [1,2,3,4,5,6,7,8,9,10,11,12,13];
                hebind = [  1,2,    3,  4,5, 6, 7, 8, 9];
                engind = [  2,3,    6,  7,8, 9,10,11,12];

                % Make paired subjects
                alleng = ENG(:,:,engind);
                allheb = HEB;
                assert(size(ENG,3),size(HEB,3))
            otherwise
                error('This case does not exist')

        end
    end

    % Calculate measures
    mval   = mean(alleng-allheb, 3);
    stdval = std(alleng-allheb, [],3);
    Cd     = zeros(128,128);
    for ii=1:128;for jj=1:128
        Cd(ii,jj)=computeCohen_d(alleng(ii,jj,:),allheb(ii,jj,:),'paired');
    end;end





    axes(ha(nr));
    imagesc(mval);axis equal;colormap(jet);colorbar;grid
    title(sprintf('%s >> Mean of diffs [EngFOV-HebFOV]',strrep(strrep(roiName,'WangAtlas_',''),'_','\_')))
    xlim([1,128]);ylim([1,128])
    xticks([1,64,128]); yticks([1,64,128])
    xticklabels([-7,0,7]); yticklabels([-7,0,7])
    caxis([-.2,.2])
    xlabel('Degs'); ylabel('Degs')

    axes(ha(nr+1*6));
    imagesc(stdval);axis equal;colormap(jet);colorbar;grid
    title(sprintf('%s >> SD of diffs [EngFOV-HebFOV]',strrep(strrep(roiName,'WangAtlas_',''),'_','\_')))
    xlim([1,128]);ylim([1,128])
    xticks([1,64,128]); yticks([1,64,128])
    xticklabels([-7,0,7]); yticklabels([-7,0,7])
    caxis([0,.45])
    xlabel('Degs'); ylabel('Degs')

    axes(ha(nr+2*6));
    imagesc(Cd);axis equal;colormap(jet);colorbar;grid
    title(sprintf("%s >> Cohen's d [EngFOV-HebFOV]",strrep(strrep(roiName,'WangAtlas_',''),'_','\_')))
    xlim([1,128]);ylim([1,128])
    xticks([1,64,128]); yticks([1,64,128])
    xticklabels([-7,0,7]); yticklabels([-7,0,7])
    caxis([-1,1])
    xlabel('Degs'); ylabel('Degs')
end
set(ha(1:12),'XTickLabel',''); set(ha,'YTickLabel','')
set(ha(1:12),'XTick',''); set(ha,'YTick','')
KK=ha(1).XLabel;
set(ha(1:12),'XLabel',KK); set(ha,'YLabel',KK)


KK.String='Degs';
set(ha(13:end),'XTickLabel',[-7,0,7]); set(ha([1,7,12]),'YTickLabel',[-7,0,7])
set(ha(13:end),'XTick',[1,64,128]); set(ha([1,7,12]),'YTick',[1,64,128])

%}

% TESTS JUST VOTRC
%{
rmroiCell_VOTRC = rmroiCell(:,6,1:2);
list_roiNames   ={'lVOTRC' };
list_dtNames     = {'Words_English','Words_Hebrew'};
list_rmNames     = {'retModel-Words_English-css-fFit.mat'
                        'retModel-Words_Hebrew-css-fFit.mat' };
% Launch the function
fname = 'Coverage_EngHeb_';  %'Fig1_'; % '' for not saving
[RF_mean, RF_individuals,empties] = figFunction_coverage_maxProfile_group(cr,list_subInds, 'flip',true, ...
                                      'rmroiCell',rmroiCell_VOTRC, ...
                                      'list_roiNames', list_roiNames, ...
                                      'list_dtNames', list_dtNames, ...
                                      'list_rmNames', list_rmNames, ...
                                      'sizedegs',7,...
                                      'minvarexp', 0.2, ...
                                      'fname', '', ...
                                      'vers',['v01_' whatFit 'fit'],...
                                      'invisible',true);
%}

% PLOT THEM FOR VOTRC, BUT LEAVING ONE SUBJECT AT A TIME
%{
% engind = [  2,3,    6,  7,8, 9,10,11,12];
ALLeng=RF_individuals{1}(:,:,engind);
ALLheb=RF_individuals{2};
% mrvNewGraphWin('EngvsHebFOVreplacement',[],true);
ha = tight_subplot(3,9,[.01 .03],[.1 .01],[.01 .01])
for kk=1:9
    % Remove 1 each time and create same plots with the remaining one
    alleng = ALLeng(:,:,[1:kk-1,kk+1:9]);
    allheb = ALLheb(:,:,[1:kk-1,kk+1:9]);

    % Calculate measures
    mval   = mean(alleng - allheb, 3);
    stdval = std(alleng  - allheb, [],3);
    Cd     = zeros(128,128);
    for ii=1:128;for jj=1:128
            Cd(ii,jj)=computeCohen_d(alleng(ii,jj,:),allheb(ii,jj,:),'paired');
    end;end

    % PLOT
    % subplot(3,9,kk) % (kk*3)-2)
    axes(ha(kk));
    imagesc(mval);axis equal;colormap(jet);colorbar;grid
    title('Mean diffs [Eng-Heb]')
    xlim([1,128]);ylim([1,128])
    % xticks([1,64,128]); yticks([1,64,128])
    % xticklabels([-7,0,7]); yticklabels([-7,0,7])
    % xlabel('Degs'); ylabel('Degs')

    % subplot(3,9,kk+9) % (kk*3)-1)
    axes(ha(kk+9));
    imagesc(stdval);axis equal;colormap(jet);colorbar;grid
    title('SD diffs [Eng-Heb]')
    xlim([1,128]);ylim([1,128])
    % xticks([1,64,128]); yticks([1,64,128])
    % xticklabels([-7,0,7]); yticklabels([-7,0,7])
    % xlabel('Degs'); ylabel('Degs')

    % subplot(3,9,kk+2*9) % kk*3)
    axes(ha(kk+2*9));
    imagesc(Cd);axis equal;colormap(jet);colorbar;grid
    title("Cohen's d")
    xlim([1,128]);ylim([1,128])
    % xticks([1,64,128]); yticks([1,64,128])
    % xticklabels([-7,0,7]); yticklabels([-7,0,7])
    % xlabel('Degs'); ylabel('Degs')
end
set(ha(1:18),'XTickLabel',''); set(ha,'YTickLabel','')
set(ha(1:18),'XTick',''); set(ha,'YTick','')
KK=ha(1).XLabel;
set(ha(1:18),'XLabel',KK); set(ha,'YLabel',KK)


KK.String='Degs';
set(ha(19:end),'XTickLabel',[-7,0,7]); set(ha([1,10,19]),'YTickLabel',[-7,0,7])
set(ha(19:end),'XTick',[1,64,128]); set(ha([1,10,19]),'YTick',[1,64,128])
set(ha(19:end),'XLabel',KK); set(ha([1,10,19]),'YLabel',KK)

%}



% PLOT THEM FOR VOTRC, DO BOOTSTRAPPING AND AVERAGE IT
% engind = [  2,3,    6,  7,8, 9,10,11,12];
% When using checkers
% old fit
%{
allind = [1,2,3,4,5,6,7,8,9,10,11,12,13];
chcind = [    1,  2,    3,4, 5, 6, 7, 8];
engind = [1,2,3,4,5,6,  7,8, 9,10,11,12];
engind = [    3  ,5,    7,8, 9,10,11,12];
%}
% new fit
allind = [1,2,3,4,5,6,7,8,9,10,11,12,13];
% R2 > 20%
chcind = [    1,  2,    3,4, 5, 6, 7   ];
engind = [1,2,3,4,5,6,  7,8, 9,10,11,12];
engind = [    3  ,5,    7,8, 9,10,12];
% R2 > 10%
chcind = [  1,2,3,4,    5,6, 7, 8, 9,10];
engind = [1,2,3,4,5,6,7,8,9,10,11,12,13];
engind = [  2,3,4,5,    8,9,10,11,12,13];


ALLeng=RF_individuals{1}(:,:,engind);
ALLheb=RF_individuals{2};

MVALS = zeros(128,128,50);
SVALS = zeros(128,128,50);
DVALS = zeros(128,128,50);

for kk=1:100
    % Remove 1 each time and create same plots with the remaining one
    assert(isequal(size(ALLeng,3),size(ALLheb,3)))
    randReplacement = datasample(1:size(ALLheb,3),size(ALLheb,3));

    alleng          = ALLeng(:,:,randReplacement);
    allheb          = ALLheb(:,:,randReplacement);

    % Calculate measures
    mval   = mean(alleng - allheb, 3);
    stdval = std(alleng  - allheb, [],3);
    Cd     = zeros(128,128);
    for ii=1:128;for jj=1:128
            Cd(ii,jj)=computeCohen_d(alleng(ii,jj,:),allheb(ii,jj,:),'paired');
    end;end

    % Accummlate it
    MVALS(:,:,kk) = mval;
    SVALS(:,:,kk) = stdval;
    DVALS(:,:,kk) = Cd;

end

% obtain means again
mval   = mean(MVALS,3);
stdval = mean(SVALS,3);
Cd     = mean(DVALS,3);

% min(Cd(:))
% max(Cd(:))
% min(Cd(:))
% mean(Cd(:))
%
% 100*sum(Cd(:)>0)/length(Cd(:))
% 100*sum(Cd(:)<=0)/length(Cd(:))


mrvNewGraphWin('CrossValEngvsHebFOV','wide',true);

caxismin = min(Cd(:));
caxismax = max(Cd(:));

subplot(1,2,1)
imagesc(Cd);axis equal;colormap(parula);colorbar;grid
title("Cohen's d: mean of 50 crossvals [Eng-CB]")
xlim([1,128]);ylim([1,128])
xticks([1,64,128]); yticks([1,64,128])
caxis([caxismin,caxismax])
xticklabels([-7,0,7]); yticklabels([-7,0,7])
xlabel('Degs'); ylabel('Degs')
% cmap = colormap(parula);
% cmap1 = cmap(1:5:128,:);
% cmap2 = cmap(129:end,:);
% vals  = linspace(0,1,50)';
% lght  = linspace(0,.3,50)';
% cmap1 = [vals,lght,lght];
% cmap2 = [lght,lght,vals];
% cmap=[flip(cmap2);cmap1];
% cmap=[cmap1;cmap2];
% cmap = cmap(129:end,:);
% colormap(cmap);
colorbar





subplot(1,2,2)
[X,Y] = meshgrid(1:128,1:128);
XX = ((X-64)/64)*7;
YY = ((Y-64)/64)*7;
YY = flipud(YY);
surf(XX,YY,Cd);
xlabel('X (degs)'); ylabel('Y (degs)')
zlabel("Cohen's d")


xlim([-7,7])
ylim([-7,7])
zlim([caxismin,caxismax])
xticks([-7,-5,-3,-1,1,3,5,7])
yticks([-7,-5,-3,-1,1,3,5,7])
xticklabels({'-7','-5','-3','-1','1','3','5','7'})
yticklabels({''})
set(gca,'FontSize',18)
caxis([caxismin,caxismax])
% colormap(cmap);
colormap(parula);











%% COVERAGE: individual plots
for subind = [1:12,14:16,18:20] % list_subInds
    subname = cr.bk.list_sub{subind}
    %% Plot the coverage figures
    % Read the coverage figure params
    % run(cr.bk.list_coverageFigure_defaults{subind});
    % Defaults
    covfig.vfc              = ff_vfcDefault;
    covfig.titleDescript    = 'FOV';
    % vfc threshold
    covfig.cothresh         = 0.2;
    covfig.vfc.cmap         = 'hot';
    covfig.vfc.addCenters   = true;
    covfig.vfc.contourPlot  = true;
    % ROIs
    covfig.list_roiNames    = {'lVOTRC'};
    covfig.list_roiNames    = {'WangAtlas_VO1_left'};
    % dt and rm names
    covfig.list_dtNames     = {'Checkers','Words','FalseFont'};
    covfig.list_rmNames     = {'retModel-Checkers-css-fFit.mat'
                               'retModel-Words-css-fFit.mat'
                               'retModel-FalseFont-css-fFit.mat'};
    % {
    % Old original fits, basically they are the same
    covfig.list_rmNames     = {'retModel-Checkers-css.mat'
                               'retModel-Words-css.mat'
                               'retModel-FalseFont-css.mat'};
    %}
    covfig.list_rmDescripts = {'Checkers', 'Words','FalseFont'};



    cr.subj.(subname).params.covfig = covfig;
    clear('covfig');
    % Plot it
    [RFcov,weight, data] = figFunction_coverage_individual(cr, subind);
    % Clean workspace of globals after each subject finishes
    mrvCleanWorkspace;
end

%% Notes
% + this is lVOTRC, do the same for V1-4,hvo1 (the ones in the paper)
% - two plots, or separate long versus short lines
% + go to white background
% - obtain numbers that show that eccc>eccw is basically noise (it will be
%      only a caption in the figure)
% + DO NOT plot any ecc diff of +- 0.5 deg
% - In the scatterplot with the light blue cones:
%    + Remove cone
%    + add +-1.5deg band
%    + below, we will only have about 3% of the data
% - compare the - and + differences, they shuold be the same for V1 v2 and
% then start changing
% - when re-running the fits, use two different HRFs (simulate that they
% will actually have size differences) and show that the effect and the
% centers will not vary
% - test: select and HRF that gives the correct size in V3 at 5deg eccen,
% and run all the analyses with this HRF
% - Use all V1, not only ventral


% test of radiality from 5 to 7 deg and 8 to 12 degs, move CB to the horizontaal and maintain the angle for words



% TODO:
% - print the line plots again  and finish figure 3
% - plot the results in IPS to check if they hold

%% DATA ANALYSIS DONE BY ROSEMARY
%{

%% (once) Upload data
% Original data in black.stanford.edu
% Create local structure that we couls upload to FW later
% - Prepare using the /utilities/prepareDataUploadFW.py script
% - Upload to FW using the fw command line utility
%      fw import folder .
%   Leonardo did it from his computer
% - TODO: add the subject and acqu metadata

%% (once) Run the analysis and upload the data table to the collection
% This example below is what we use for DWI analysis with RTP. Now we need
% something similar for PRF-s: run them in the server, and download only the
% relevant data for our analyses and plots (and next we could obtain the
% analyses online as well.

%   Usually I run this in Google server and then continue locally
%{
% Run the analyses
serverName     = 'stanfordlabs';
collectionName = 'HCP_Depression';
gearName       = 'afq-pipeline';
gearVersion    = '3.0.7';
% Before launching check the script for the correct analysis config params
dr_fwLaunchJobs(serverName, collectionName, gearName, gearVersion)

% Create the datatable (only fa for now)
measurement    = 'fa';
% get all the analysis in the collection
JL = dr_fwCheckJobs(serverName, collectionName);
% filter the analyses
state         = 'complete';
dateFrom      = '04-Feb-2019 00:00:00';
labelContains = 'v3.0.7:';
t  = JL(JL.state==state & JL.gearName==gearName & ...
       JL.gearVersion==gearVersion & JL.JobCreated>dateFrom & ...
       contains(string(JL.label), labelContains),:);
% read data of interest and create table
dt = dr_fwReadDtFromAnalysisTable(serverName, t, measurement);
% It takes a lot of time, save it locally...
localfname    = fullfile(stRootPath,'local','tmp', ...
                      sprintf('AllV02_HCP_Depression_%s.mat',measurement));
save(localfname, 'dt')
%  ...and upload it to the collection
st   = scitran(serverName); st.verify;
cc   = st.search('collection','collection label exact',collectionName);
stts = st.fileUpload(localfname, cc{1}.collection.id, 'collection');
% Check that the data is there
[~,fname,ext] = fileparts(localfname);
data          = load(st.fw.downloadFileFromCollection(cc{1}.collection.id,...
                                                  [fname ext],localfname));
%}

%% Download the data from the collection for analysis
%{
% Every time we want to re-run the data, we can download it from the server.

% Download or clone the repository
%       !git clone https://github.com/garikoitz/paper-HCPDEPRESSION.git
%
clear all; close all; clc;
% Add the root of the repository to the Matlab path (or run this code):

    % cd('<path-to-your-code>/paper-reproducibility')
    cd('~/soft/paper-HCPDEPRESSION')
    rootDir = pwd;
    addpath(genpath(rootDir));

% Specify a path to save the output figures.
paperPath = '~/gDrive/STANFORD/PROJECTS/2019 Depression and WM (Leonardo-Gari-Brian)';
saveItHere = string(fullfile(paperPath, 'VERSION_01/figures/sources'));


% Read the data
% Check if there is a local cache, otherwise download it from FW
DataVersion    = '01';
collectionName = 'HCP_Depression';
measure        = 'fa';

fname          = sprintf('AllV%s_%s_%s.mat',DataVersion, collectionName, measure);
localfname = fullfile(paperReprPath,'local',fname);
if exist(localfname,'file')
    data = load(localfname);
else  % Download it from the Flywheel collection attachment
    serverName     = 'stanfordlabs';
    st  = scitran(serverName);
    cc  = st.search('collection','collection label contains',collectionName);
    data=load(st.fw.downloadFileFromCollection(cc{1}.collection.id,fname,localfname));
end

% This script can be run directly using the Run button or step by step.
%}

%% Data preparation
if (0)
% Do not do it, use what we have already

% apply canonical x form -- for every nifti
niftiApplyCannonicalXform

% acpc align the anatomical
mrAnatAverageAcpcNifti

% run freesurfer
eval(['! recon-all -i ' pathT1 ' -subjid ' dirNameFreesurfer ' -all'])

% ribbon from freesurfer into class file -- t1_class.nii.gz
fs_ribbon2itk(inputRibbonFile, outputClassNii, [], pathT1, [])

% TODO: see if we can substitute the previous steps using fmriprep
end

%% Initialization
if (0)
% Do not do it, use what we have already

% Initialize. This will create the mrSESSION in the root folder.
cr_mrInit(cr, opt);
end

%% Allignment and dataType creation
if (0)
% Do not do it, use what we have already

% align the inplane to anatomical
s_alignInplaneToAnatomical;

% This is the message when I closed it:
%    The alignvolumedata GUI has been closed. in case you forgot to export the transformation, here it is:
%    tr = maketransformation([0 0 0],[1 2 3],[95.9343856857804 106.096638132966 92.1256300034666],[1 2 3],[-88.7924000052859 -2.12048422951174 98.6511462809804],[192 192 62],[208.000007629395 208.000007629395 62],[0.998614818358445 -1.00108310465201 1.9968769575349],[0 0 0],[0.000745334952393827 -0.000334517493656755 0.000639298200901277],[0 0 0]);
% make the transformation into a 4x4 matrix
%    rx.volVoxelSize = [2 2 2];
%    T = transformationtomatrix(tr,0,rx.volVoxelSize);
% vw = initHiddenInplane; mrGlobals;
% mrSESSION.alignment = T;
% saveSession;




% TODO??
% specify segmentation file, go to gray view to run the prfs

% create a new dataTYPE which is the average of the 4 runs (note, still in INPLANE)
% this script will also xform the data into the gray

% the names of the dataTYPES we want to create
dtsToCreate = {
    'WordFalse1';     % 1 % the 1st checker and the 2nd word
    'WordFalse2';     % 2 % the 2nd checker and the 1st word
    };
opt.dtsToCreate = {
    'Checkers1'         % 1
    'Words_English1'    % 2
    'Words_Hebrew1'     % 3
    'Words_Hebrew2'     % 4
    };
%
% The datatype the scan belongs to. For example, a 1 means that the first
% scan is in the first dataTYPE specified in dtsToCreate
% use 0 if the scan is not used in the creation of a new dt
opt.dtAssignments = [
    0;
    0;
    1;
    2;
    3;
    4;
    ];

% make the new tseries from the most processed time series
opt.dtToAverage = 'MotionComp_RefScan1';
% Run it
cr = crNewDataTypes(cr, opt);
end

%% Stimuli preparation
if (0)
% Do not do it, use what we have already

% make a Stimuli folder in the same place as the mrSESSION.mat
% for localizer GLM analyses, make Stimuli/Parfiles
% 2 things go into Stimuli
% params file -- mrVista writes this to desktop
% image matrix -- (part of the params file)
end

%% Analysis of the data
if (0)
% The code here will go to a Docker container. Make the figures in the paper reproducible
% The code below will come as a combination from pmVistasoft.m that I did
% adapting a script from Jon Winawer, and s_prfRun.m Rosemary Le's script.

cr.dirs.prfRun_params()



% subjects we want to do this for
opt.list_subInds = [31];
% dataTYPE name. Can run for mutiple datatypes
opt.list_rmName = {'Words_Hebrew1','Words_English1'};
% roi name. assumes in shared anatomy directory (change this to be self contained)
% if we want to run on the whole brain, assign this the empty string ''
% assign this to be a string in a cell otherwise {'LV1_rl'}
opt.list_rois = {'lVOTRC'};
% prf model. Specify in a cell. Options:
% {'one oval gaussian' | 'onegaussian' | 'css'}
% Note: if we want to specify multiple models, change the naming
% convention. See outFileName
opt.prfModel = {'css'};
% search type.
% 1 = grid search only ("coarse"),
% 2 = minimization search only ("fine"),
% 3 = grid followed by minimization search [default]
%   note, there is another option which is to find the hrf as well
opt.wSearch = 3;
% radius of circle retinotopy in visual angle degrees
opt.p.stimSize = 7; %%% Is this true for the Hebrew subject we just selected?
% define things common to all datatypes
% name of params file
opt.p.paramsFile_Knk        = 'Stimuli/params_knkfull_multibar_blank.mat';  % Words and FalseFont
opt.p.paramsFile_Checkers   = 'Stimuli/params_checkers.mat';                % Checkers
% image file
opt.p.imFile_Knk            = 'Stimuli/images_knk_fliplr.mat';              % Words and FalseFont
opt.p.imFile_Checkers       = 'Stimuli/images_8barswithblank_fliplr.mat';   % Checkers
% params common to all dts
opt.params.stimSize         = opt.p.stimSize;
opt.params.fliprotate       = [0 0 0];
opt.params.stimType         = 'StimFromScan';
opt.params.stimWidth        = 90;
opt.params.stimStart        = 0;
opt.params.stimDir          = 0;
opt.params.nCycles          = 1;
opt.params.nStimOnOff       = 0;
opt.params.nUniqueRep       = 1;
opt.params.nDCT             = 1;
opt.params.hrfType          = 'two gammas (SPM style)';
opt.params.hrfParams        = {[1.6800 3 2.0500] [5.4000 5.2000 10.8000 7.3500 0.3500]};
opt.params.imfilter         = 'binary';
opt.params.jitterFile       = 'Stimuli/none';

% Run it

% We downloaded the mrSESSION from black, edit it before using it locally
results = cr_prfRun(cr, opt);
end

%% Generate the individual figure
if (0)
% {
% For a full list of options, see edit IndivFigDefaults
% title description
opt.titleDescript    = 'FOV';
% vfc threshold
opt.vfc              = ff_vfcDefault_Hebrew;
opt.vfc.cmap         = 'hot';
opt.vfc.addCenters   = true;
opt.vfc.contourPlot  = true;
% subjects
opt.list_subInds     =  [31];
% session
% opt.list_path        = cr.bk.list_sessionRet;
% opt.list_path        = list_sessionRet;
% ROIs
opt.list_roiNames    = {'lVOTRC'};
% dt and rm names
opt.list_dtNames     = {'Words_Hebrew','Words_English'};
opt.list_rmNames     = {'retModel-Words_Hebrew-css.mat'
                        'retModel-Words_English-css.mat'};
opt.list_rmDescripts = {'Words_Hebrew', 'Words_English'};
figFunction_coverage_individual(cr, opt);
%}
% end

%}

%% Line plots: from checquerboard to word
% Combine code in figScript_centers_* to create the plots Brian nd I
% discussed
%   - left VOTC, line plot going from word to cb
%   - black if both word and cb are in the same quadrant (); else gray
%   - continuous line if word is inside fovea (0-x (x=2deg)), and CB is out; else dashed
%   - count percentage of continuous black lines over rest:
%      - bin by size (small sizes very unreliable)
%      -

% SELECT SUBJECTS AND MODELS
%{
onlyStanford  = [1:20];  % Why not the rest? Ask MBS/RL
onlyHebrew    = [31:36 38:44];
list_subInds  = [onlyStanford onlyHebrew];
list_subInds  = [1:12];
list_subInds = [1:12,14:16,18:20];

list_path     = cr.bk.list_sessionRet;
list_roiNames = {'WangAtlas_V1v_left';
                 'WangAtlas_V2v_left';
                 'WangAtlas_V3v_left';
                 'WangAtlas_hV4_left';
                 'WangAtlas_VO1_left';
                 'lVOTRC';
                 'WangAtlas_IPS0';
                 'WangAtlas_IPS1'};
% list_roiNames = {'LV1_rl'
%                  'LV2v_rl'
%                  'LV3v_rl'
%                  'LhV4_rl'
%                  'LVO1_rl'
%                  'lVOTRC' };
list_dtNames = {'Words'...
                'Checkers'...
                'Words_English'...
                'Words_Hebrew'...
                'FalseFont'};
% ret model names
list_rmNames = {'retModel-Words-css-fFit.mat'...
                'retModel-Checkers-css-fFit.mat'...
                'retModel-Words_English-css.mat'...
                'retModel-Words_Hebrew-css.mat'...
                'retModel-FalseFont-css.mat'};
list_rmDescripts = {'Words'...  % Words (large bars)
                    'Checkers'...
                    'Words_English'...
                    'Words_Hebrew'... % Words (smalls bars)
                    'FalseFont'};

%% Get the cell of rms so that we can threshold
if 0
    rmroiCell = ff_rmroiCell(cr,list_subInds, list_roiNames, list_dtNames, ...
            list_rmNames, 'list_path', list_path);
    % SAVE THIS TO WORK LOCALLY
    mkdir(fullfile(sdRP,'DATA'))
    % 20 subs, words and CBs
    % save(fullfile(sdRP,'DATA','rmroicell_1to20.mat'),'rmroiCell');
    % ALL subs, words and CBs and FF and Heb
    save(fullfile(sdRP,'local','rmroicell_allandall.mat'),'rmroiCell');
end
if 0
    % load(fullfile(sdRP,'DATA','rmroicell_1to20.mat'));
    load(fullfile(sdRP,'local','rmroicell_allandall.mat'));
    % Only the 20
    list_subInds  = [1:20];
    rmroiCell = rmroiCell(list_subInds,(1:6),[1,2]);
    % Only the ones with FF
    % list_subInds  = [1:12];
    % rmroiCell = rmroiCell(list_subInds,(1:6),[1,5]);
    % Only the Hebrew
    % rmroiCell = rmroiCell((21:end),(1:6),[1,5]);
end

%%  Threshold and get identical voxels for each subject
% values to threshold the RM struct by
% vfc = ff_vfcDefault_Hebrew;
% The defaults are different for the two projects, sofor now use the most
% restrictive one
% vfc threshold
vfc.prf_size         = true;
vfc.fieldRange       = 15; % 7;
vfc.method           = 'max';
vfc.newfig           = true;
vfc.nboot            = 50;
vfc.normalizeRange   = true;
vfc.smoothSigma      = false;
% Thresholds
vfc.sigmaEffthresh   = [.2 15];  % [0  7]; % sigma effect (sigmaMajor/sqrt(exponent))
vfc.sigmaMajthresh   = [0.5 8];   % [0 14]; % sigma major (before the exponent)
vfc.cothresh         = 0.20; % Variance explained, 20%
vfc.cothreshceil     = 1; % 0.2 looking for noise. 1 for normal thresholding. don't get voxels higher than this
vfc.threshByCoh      = true;
vfc.eccthresh        = [0 vfc.fieldRange];
vfc.quadthresh       = [1 4];  % If center not on these quads, remove
% Fig control
vfc.nSamples         = 128;
vfc.meanThresh       = 0;
vfc.weight           = 'fixed';
vfc.weightBeta       = false;
vfc.cmap             = 'jet';
vfc.clipn            = 'fixed';
vfc.addCenters       = false;
vfc.verbose          = prefsVerboseCheck;
vfc.dualVEthresh     = 0;
vfc.ellipsePlot      = false;
vfc.ellipseLevel     = 0.5;
vfc.ellipseColor     = [1 0 0];
vfc.contourPlot      = true;
vfc.contourLevel     = 0.5;
vfc.contourColor     = [0 0 0];
vfc.tickLabel        = false;
vfc.contourBootstrap = false;
vfc.gridColor        = [.6 .6 .6];
vfc.backgroundColor  = [1 1 1];
vfc.cmapRange        = [0 pi]; % the range over which color bar
vfc.cmapValues       = flipud(jetCmap(0,128)); % the colorbar values
vfc.alphaValue       = ''; % 0.5;
vfc.alphaValueDot    = ''; % 0.8;
vfc.lineWidth        = 1.5; % thicker lines for transparent. 1 works well for opaque
% by definition, when eccentricity does not shift a lot, theta will not
% shift either. We can look at theta shifts in the voxels whose
% eccentricity have shifted a lot
% thetaShiftByEccThresh = true;
% eccThresh = 3;
% Update the cr structure that I created that go into the new functions
cr.defaults.covfig.vfc = vfc;
% INITIALIZE SOME THINGS
numRois = length(list_roiNames);
numSubs = length(list_subInds);
% cell for linearizing the data (a vector for each ROI)
L_data = cell(1, numRois);
X_rm1 = cell(1, numRois);
Y_rm1 = cell(1, numRois);
X_rm2 = cell(1, numRois);
Y_rm2 = cell(1, numRois);
Ecc_rm1 = cell(1, numRois);
Ecc_rm2 = cell(1, numRois);
rmDescript1 = list_rmDescripts{1};
rmDescript2 = list_rmDescripts{2};
% In comparing ret models, the collection of voxels may not be the same
% because of the thresholding. In this cell we redefine the rmroi
rmroiCellSameVox = cell(size(rmroiCell));
for jj = 1:numRois
    for ii = 1:numSubs
        % get identical voxels for each subject's roi over all ret models
        D = rmroiCell(ii,jj,:);
        % GLU EDIT function: remove voxels from the oppossite hemifield
        rmroiCellSameVox(ii,jj,:) = ff_rmroiGetSameVoxels(D, cr.defaults.covfig.vfc);
    end
end

% Linearize the data
% Take the difference between 2 rms.
% Also store the x and y data
for jj = 1:numRois
    % initializing the difference of the centers' thetas
    ldata = [];

    % intializing the location of the centers
    xdata_rm1 = [];
    ydata_rm1 = [];
    xdata_rm2 = [];
    ydata_rm2 = [];

    % initializing eccentrcity
    ecc_rm1   = [];
    ecc_rm2   = [];

    % initializing angle
    ph_rm1    = [];
    ph_rm2    = [];

    % initializing size
    sm_rm1    = [];
    sm_rm2    = [];

    for ii = 1:numSubs
        rmroi1 = rmroiCellSameVox{ii,jj,1};
        rmroi2 = rmroiCellSameVox{ii,jj,2};

        % some subjects don't have
        if ~isempty(rmroi1) & ~isempty(rmroi2)
            data1 = rmroi1.ph;
            data2 = rmroi2.ph;

            % the difference between centers' thetas.
            % this will determine the color of the line
            % we take absolute value because we are interested in the magnitude
            % of the rotation and not the direction
            fieldDiffOver = abs(data2 - data1);

            % Note that the difference will range between 0 and 2pi.
            % We want to constrain values to be between and pi (again not
            % interested in the direction of the rotation but the magnitude of it)
            % For values greater than pi, subtract it from 2pi
            fieldDiff = ff_polarAngleBetween0AndPi(fieldDiffOver);

            ldata = [ldata fieldDiff];

            % the location of the pRF centers
            xdata_rm1 = [xdata_rm1 rmroi1.x0];
            ydata_rm1 = [ydata_rm1 rmroi1.y0];

            xdata_rm2 = [xdata_rm2 rmroi2.x0];
            ydata_rm2 = [ydata_rm2 rmroi2.y0];

            ecc_rm1   = [ecc_rm1 rmroi1.ecc];
            ecc_rm2   = [ecc_rm2 rmroi2.ecc];

            ph_rm1    = [ph_rm1 rmroi1.ph];
            ph_rm2    = [ph_rm2 rmroi2.ph];

            sm_rm1    = [sm_rm1 rmroi1.sigma1];
            sm_rm2    = [sm_rm2 rmroi2.sigma1];
        end
    end
    L_data{jj} = ldata;

    X_rm1{jj}   = xdata_rm1;
    Y_rm1{jj}   = ydata_rm1;

    X_rm2{jj}   = xdata_rm2;
    Y_rm2{jj}   = ydata_rm2;

    Ecc_rm1{jj} = ecc_rm1;
    Ecc_rm2{jj} = ecc_rm2;

    Ph_rm1{jj}  = ph_rm1;
    Ph_rm2{jj}  = ph_rm2;

    Sm_rm1{jj}  = sm_rm1;
    Sm_rm2{jj}  = sm_rm2;

end
% Get a colormap according to the linearized data in L_data
for jj = 1:numRois
    ldata = L_data{jj};
    cdata = ff_colormapForValues(ldata, cr.defaults.covfig.vfc.cmapValues, ...
                                        cr.defaults.covfig.vfc.cmapRange);
    C_data{jj} = cdata;
end

% colormap for histogram
% cmapValuesHist = colormap('pink');
% cmapValuesHist_tem = colormap('hot');
% cmapValuesHist = cmapValuesHist_tem(2:55, :);
colormap(zeros(64,3)); % matlab has funky behavior where the size of this influences the size of all future colorbars...
cmapValuesHist = colormap('pink');
close;

%% PLOT IT
for fov      = [ 1.5]
xx = mrvNewGraphWin('LineRadiality and Scatterplot','wide');
set(xx,'Position',[0.005 0.062 .95 .55 ]);
for jj = 1:numRois
    % data
    ldata = L_data{jj};

    X1 = X_rm1{jj};
    Y1 = Y_rm1{jj};

    X2 = X_rm2{jj};
    Y2 = Y_rm2{jj};

    C = C_data{jj};

    ecc_rm1 = Ecc_rm1{jj};
    ecc_rm2 = Ecc_rm2{jj};

    ph_rm1 = Ph_rm1{jj};
    ph_rm2 = Ph_rm2{jj};

    sm_rm1 = Sm_rm1{jj};
    sm_rm2 = Sm_rm2{jj};

    roiName = list_roiNames{jj};


    % plot on polar map
    %   + black if both word and cb are in the same quadrant (); else gray
    %   + continuous line if word is inside fovea (0-x (x=2deg)), and CB is out; else dashed
    %   + count percentage of continuous black lines over rest:
    %      + bin by size (small sizes very unreliable)

    % Filter by size, always, +-inf for no filter
    % Moved this threshold above, now I can't count, make reports in the f
    % sizeMIN  = .5;
    % sizeMAX  = 8;

    % sIndw = (sm_rm1 <= sizeMAX) & (sm_rm1 >= sizeMIN);
    % sIndc = (sm_rm2 <= sizeMAX) & (sm_rm2 >= sizeMIN);
    % sInd  = (sIndw & sIndc);

    % ldata = ldata(sInd);
    % X1 = X1(sInd);
%     Y1 = Y1(sInd);
%     X2 = X2(sInd);
%     Y2 = Y2(sInd);
%     C = C(sInd);
%     ecc_rm1 = ecc_rm1(sInd);
%     ecc_rm2 = ecc_rm2(sInd);
%     ph_rm1  = ph_rm1(sInd);
%     ph_rm2  = ph_rm2(sInd);
%     sm_rm1  = sm_rm1(sInd);
%     sm_rm2  = sm_rm2(sInd);
%     % histogram(sm_rm1);hold on;histogram(sm_rm2);legend()
    % Obtain angle and ecc again, just in case
    % [PW,ECCW] = cart2pol(X1,Y1);
    % [PC,ECCC] = cart2pol(X2,Y2);

   if 0
        fovs   = [0,0.5,1,1.5, 2,3,4,5];
        for ff = 1:length(fovs)
            fov = fovs(ff);
            counter = 0;
            for pp = 1:length(X1)
                pw   = rad2deg(PW(pp));
                pc   = rad2deg(PC(pp));
                qw   = floor(pw/90)+1;
                qc   = floor(pc/90)+1;
                eccw = ECCW(pp);
                eccc = ECCC(pp);

                % if qw==qc && eccw < fov && eccc > fov; counter=counter+1; end
                % if qw==qc && (eccc - eccw) > fov; counter=counter+1; end
                % if qw==qc && (eccc > eccw); counter=counter+1; end
                if (-eccc + eccw) > fov; counter=counter+1; end
                % if eccw < fov && eccc > fov; counter=counter+1; end
                % if qw==qc && ((eccc - eccw) > fov) ; counter=counter+1; end
            end
            if ff==1
                fprintf('\n\n--------------------------------------------------\n');
                % fprintf('SizeMIN: %g, sizeMAX: %g', sizeMIN, sizeMAX);
                % fprintf('  (orig: %g, filtered: %g)', length(sInd), sum(sInd));
                fprintf('\n---------------------------------------------------\n');
            end
            fprintf('Same quadrant and word inside %g deg and cb out, %i out of %i, %0.2f%%\n',...
                    fov,counter,length(X1),counter*100/length(X1))
        end
   end
    testing = false;
    radialityCone=true;
   if 1
        % initialize polar plot
        % Limit plot to visual field circle
        subplot(1,numRois,jj);
        axis([-cr.defaults.covfig.vfc.fieldRange cr.defaults.covfig.vfc.fieldRange ...
              -cr.defaults.covfig.vfc.fieldRange cr.defaults.covfig.vfc.fieldRange])

        % polar plot
        ff_polarPlot(cr.defaults.covfig.vfc);
        hold on;

        % colorbar
        % c = colorbar;
        % colormap(cr.defaults.covfig.vfc.cmapValues)
        % set(c, 'Color', [1 1 1])
        % caxis(cr.defaults.covfig.vfc.cmapRange)
        counter  = 0;

        maxang15 = [];
        minang15 = [];
        maxang5  = [];
        minang5  = [];
        for pp = 1:length(X1)
            % pp = 555;
            % lineColor = C(pp,:);
            % pw   = rad2deg(PW(pp));
            % pc   = rad2deg(PC(pp));
            % qw   = floor(pw/90)+1;
            % qc   = floor(pc/90)+1;
            % eccw = ECCW(pp);
            % eccc = ECCC(pp);

            pw   = ph_rm1(pp);
            pc   = ph_rm2(pp);
            eccw = ecc_rm1(pp);
            eccc = ecc_rm2(pp);



            % if qw==qc ; lineColor = [1 1 1];
            % else; lineColor = [.5 .5 .5]; end

            % if eccw < 5 && eccc > 5 ; lineStyle = '-';
            % else; lineStyle = ':'; end

            % if qw==qc && eccw < 5 && eccc > 5; counter=counter+1; end
            if (abs(eccc - eccw)  > fov)
                if  ((eccc - eccw) > 0) % qw==qc &&
                    counter=counter+1;
                    lineColor = [.1 .1 .1];
                    lineStyle = ['-'];
                else
                    lineColor = [.5 .5 .5];
                    lineStyle = [':'];
                end

                if radialityCone
                    % WORDS
                    x1 = X1(pp);
                    y1 = Y1(pp);
                    % CB
                    x2 = X2(pp);
                    y2 = Y2(pp);
                    % We know the polar angle and the ecc, rotate & translate
                    % Rotate angle
                    pcn    = 0;
                    pwn    = pw - pc;
                    % Convert to cartesian
                    [x1n, y1n] = pol2cart(pwn, eccw);
                    [x2n, y2n] = pol2cart(pcn, eccc);
                    % Do the translation only in the x axis
                    if eccc > 5
                        x2nt  = 15;
                        x2dif = x2nt - x2n;
                        x1nt  = x1n + x2dif;
                        % Plot with the same y
                        line([x1nt x2nt], [y1n, y2n], ...
                            'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth',1);
                    else
                         x2nt  = 5;
                         x2dif = x2nt - x2n;
                         x1nt  = x1n + x2dif;
%                         % Plot with the same y
%                         line([x1nt x2nt], [y1n, y2n], ...
%                             'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth',1);
                    end

                     % Calculate angles
                     angle = atand(y1n/(x2nt-x1nt));
                    if  ((eccc - eccw) > fov)
                        if eccc > 5
                            if angle > 0
                                maxang15 = [maxang15 angle];
                            end
                            if angle < 0
                                minang15 = [minang15 angle];
                            end
                        else
                            if angle > 0
                                maxang5 = [maxang5 angle];
                            end
                            if angle < 0
                                minang5 = [minang5 angle];
                            end
                        end
                    end

                    if testing
                        line([x1 x2], [y1, y2], ...
                            'Color', 'r', 'LineStyle', lineStyle, ...
                             'LineWidth',1);
                         % Assertions that the calculations are right
                        d      = sqrt((x2-x1)^2 +  (y2-y1)^2)
                        dn     = sqrt((x2nt-x1nt)^2 +  (y2n-y1n)^2)
                    end

                else
                    plot([X1(pp) X2(pp)], [Y1(pp), Y2(pp)], ...
                        'Color', lineColor, 'LineStyle', lineStyle, ...
                         'LineWidth',1);
                end
            end
        end

        % titleName = {
            % ['pRF Center radiality in ' roiName]
            % ['white(' sprintf('%0.2f%% ',counter*100/length(X1)) ...
            %  ')=same quadrant, ecc: (' rmDescript2 ' - ' rmDescript1 ') > ' num2str(fov)]
            % ['sMIN:' num2str(sizeMIN) ',sMAX:' num2str(sizeMAX) ...
            % ' (pRFs ' num2str(sum(sInd)) ' of '  num2str(length(sInd)) ')']
        %     };
        titleName = {
            [strrep(ff_stringRemove(roiName, 'WangAtlas_'),'_','\_')];
            % ['slope: ' num2str(meanSlope)];
            % ['ci: ' num2str(ci')];
            % [num2str(numVoxels) ' voxels']
            };
        ttext = strrep(titleName{1},'\_left','');

        % title(titleName, 'fontweight', 'bold', 'color', [.1 .1 .1], 'fontsize', 14);
        text(-9, 11, ttext,'Color','k','FontWeight','Bold','FontSize',16)

        if ~isempty(maxang15)
            text(-9, 3,sprintf('+%1.1f°',max(maxang15)),'Color','k','FontSize',10)
            text(-9, 1,sprintf('+%1.1f°',median(maxang15)),'Color','r','FontSize',10)
            line([15+5*cosd(180-median(maxang15)),15],[5*sind(median(maxang15)),0],'Color','r')
        end
        if ~isempty(minang15)
            text(-9, -1,sprintf('%1.1f°',median(minang15)),'Color','r','FontSize',10)
            text(-9, -3,sprintf('%1.1f°',min(minang15)),'Color','k','FontSize',10)
            line([15+5*cosd(180-median(minang15)),15],[5*sind(median(minang15)),0],'Color','r')
        end

%         if ~isempty(maxang5)
%             text(-9, 3,sprintf('+%1.1f°',max(maxang5)),'Color','k','FontSize',10)
%             text(-9, 1,sprintf('+%1.1f°',median(maxang5)),'Color','r','FontSize',10)
%             line([5+5*cosd(180-median(maxang5)),5],[5*sind(median(maxang5)),0],'Color','r')
%         end
%         if ~isempty(minang5)
%             text(-9, -1,sprintf('%1.1f°',median(minang5)),'Color','r','FontSize',10)
%             text(-9, -3,sprintf('%1.1f°',min(minang5)),'Color','k','FontSize',10)
%             line([5+5*cosd(180-median(minang5)),5],[5*sind(median(minang5)),0],'Color','r')
%         end
        xlim([-15, 18])
        % titlefile = strrep(titleName{1},' ','_');
        % saveas(gcf, fullfile(sdRP,'local','png',[titlefile '.png']), 'png')
%         subplot(2,numRois,numRois+jj);
%         polarhistogram(deg2rad([maxang5 minang5]),100, 'Normalization','probability',...
%             'EdgeAlpha',1, 'EdgeColor','k','FaceAlpha',1, 'FaceColor','k');
%         hold on;
%         polarhistogram(deg2rad([maxang15 minang15]),100, 'Normalization','probability',...
%             'EdgeAlpha',1, 'EdgeColor','b','FaceAlpha',.4, 'FaceColor','b');
%         legend({'< 5°','> 5°'})


   end
end
saveas(gcf, fullfile(sdRP,'local','png',['Radiality_Plots_2x' num2str(fov) '.png']), 'png')
saveas(gcf, fullfile(sdRP,'local','svg',['Radiality_Plots_2x' num2str(fov) '.svg']), 'svg')
close all
end

%% SCATTERPLOTS


% whether looking at a subject by subject basis
subIndividually = false;

% field to plot. Ex:
% 'co': variance explained
% 'ecc': eccentricity
% 'sigma': effective size
% 'sigma1': sigma major
% 'numvoxels' for number of voxels in roi
% fieldToPlotDescript is for axis labels and plot titles
%     'sigma'       : effective sigma
%     'sigma1'      : sigma major
%     'ecc'         : eccentricity
%     'co'          : variance explained
%     'exponent'    : exponent
%     'betaScale'   : how much to scale the predicted tseries by
%     'meanMax'     : mean of the top 8 values
%     'meanPeaks'   : mean of the outputs of matlab's meanPeaks
list_fieldNames  = {
    'co'
    'ecc'
%     'sigma'
%     'ph'
    };

list_fieldDescripts = {
    'variance explained'
    'eccentricity'
%     'sigma'
%     'polar angle'
    };

% which plots do we want? lots we can make ...
plot_fit = false; % plotting the across-subject bootstrapped line w/ CIs

% transparency of the plots
alphaValue = 0.4;



% location of the colorbar
% default: 'eastoutside'
% 'southoutside':
cbarLocation = 'eastoutside';

% end modification section

numSubs = length(list_subInds);
numRois = length(list_roiNames);
numRms = length(list_rmNames);

% number of fields
numFields = length(list_fieldNames);

% rm descriptions
rm1Descript = list_rmDescripts{1};
rm2Descript = list_rmDescripts{2};

% initialize structs / matrices for mixed effects
subjectLines = cell(numSubs, numRois, numFields); % because pairwise

% initialize struct for calculating the percentage of voxels above the
% identity line
percentAbovePooled = zeros(numRois, numFields);
percentAboveSubs   = zeros(numSubs, numRois, numFields);
A = cell(numFields*numRois, 5);

% get the cell of rms so that we can threshold
% rmroiCell = ff_rmroiCell(cr, list_subInds, list_roiNames, list_dtNames, list_rmNames, ...
%     'list_path', list_path);

% Threshold and get identical voxels for each subject
% In comparing ret models, the collection of voxels may not be the same
% because of the thresholding. In this cell we redefine the rmroi
% rmroiCellSameVox = cell(size(rmroiCell));

% for jj = 1:numRois
%     for ii = 1:numSubs
        % get identical voxels for each subject's roi over all ret models
%         D = rmroiCell(ii,jj,:);
%         rmroiCellSameVox(ii,jj,:) = ff_rmroiGetSameVoxels(D, vfc);
%     end
% end

% close all;


    % subplot(2,numRois,numRois+jj);


for ff = 1:numFields

    % field-specific properties
    fieldName = list_fieldNames{ff};
    fieldNameDescript = list_fieldDescripts{ff};

    if strcmp(fieldName, 'sigma1')
        maxValue = cr.defaults.covfig.vfc.sigmaMajthresh(2);
    elseif strcmp(fieldName, 'sigma')
        maxValue = cr.defaults.covfig.vfc.sigmaEffthresh(2);
    elseif strcmp(fieldName, 'ecc')
        maxValue = cr.defaults.covfig.vfc.eccthresh(2);
        fov = 1.5; % width of the band
        nrows = 2; ncols = 3;
        position = [0.005 0.062 .9 .7 ];
    elseif strcmp(fieldName, 'co')
        maxValue = 1;
        fov = 0.2; % width of the band
        nrows = 1; ncols = 6;
        position = [0.005 0.062 .9 .5 ];
    elseif strcmp(fieldName, 'exponent')
        maxValue = 2;
    elseif strcmp(fieldName, 'meanMax')
        maxValue = 20;
    elseif strcmp(fieldName, 'meanPeaks')
        maxValue = 10;
    elseif strcmp(fieldName, 'betaScale')
        maxValue = 5;
    elseif strcmp(fieldName, 'x0') || strcmp(fieldName, 'y0')
        maxValue = cr.defaults.covfig.vfc.fieldRange;
    elseif strcmp(fieldName, 'ph')
        maxValue = 2*pi;
    else
        error('Define the maxValue so we can normalize and fit the beta distribution.');
    end

    xx = mrvNewGraphWin('Scatterplots');
    set(xx,'Position',position);
    for jj = 1:numRois
        roiName = list_roiNames{jj};
        subplot(nrows,ncols,jj);



        axisLims = [0 maxValue];
        BarData1 = [];
        BarData2 = [];

        for ii = 1:numSubs

            subInd = list_subInds(ii);

            % rmRois for different ret models
            rmroi1 = rmroiCellSameVox{ii,jj,1};
            rmroi2 = rmroiCellSameVox{ii,jj,2};

             if ~isempty(rmroi1)

                % get the data
                x1 = eval(['rmroi1.' fieldName]);
                x2 = eval(['rmroi2.' fieldName]);

                % shift so that the smallest value is 0
                if strcmp(fieldName, 'ph')
                    x1 = x1 + pi;
                    x2 = x2 + pi;
                end

                %% For the mixed effects
                % fit a line for each subject
                % p = polyfit(x1, x2, 1);
                % subjectLines{ii,jj,ff} = p;
                % [B, BINT] = regress(x2', x1');
                % b.slope     = B;
                % b.ci95      = BINT;

                % p = polyfit(x1, x2, 1);
                % b.pintercept = p(1);
                % b.pslope = p(2);

                % subjectLines{ii,jj,ff} = b;

                %% the percentage of voxels above the identityLine
                perAbove = sum(x2 > x1) / length(x2);
                percentAboveSubs(ii,jj,ff) = perAbove;

                %% concatenate
                BarData1 = [BarData1, x1];
                BarData2 = [BarData2, x2];

            end
        end % end loop over subjects

        %% mixed effects: fit a line to individual subjects
        % slopes = nan(1, numSubs);
        % slopesp = nan(1,numSubs);
        % interceptsp = nan(1,numSubs);
        % percents = percentAboveSubs(:,jj,ff);

        % for ii = 1:numSubs
        %     b = subjectLines{ii,jj,ff};
        %     if ~isempty(b)
        %         slopes(ii) = b.slope;
        %     end
        % end

        % the calculating. nan will cause bootci to error
        % slopes(isnan(slopes)) = [];
        % percents(isnan(percents)) = [];

        % table things.
        % (1)roiName (2)fieldName (3) ciLow (4)ciHigh (5)mean
        % tind = (jj-1)*numFields + ff;


        % if numSubs > 1
            % numbs = 1000;
            % [ci, bootstat] = bootci(numbs, @mean, slopes);
            % meanSlope = mean(bootstat);

            % % table things for percent above identityLine
            % [ciPer, bootstatPer] = bootci(numbs, @mean, percents);
            % A{tind,1} = ff_stringRemove(roiName, '_rl');
            % A{tind,2} = fieldName;
            % A{tind,3} = ciPer(1);
            % A{tind,4} = ciPer(2);
            % A{tind,5} = mean(bootstatPer);

        % else
        %     meanSlope = nan;
        %     ci = nan;
        %   end

        %% calculations related to the percentage above the identityLine
        % percentabove = sum(BarData2 > BarData1) / length(BarData1);
        % percentAbovePooled(jj,ff) = percentabove;

        % properties related to both types of scatter plots
        % coloring by number of voxels or percentage of voxels
        npoints = 100;

        % 3d histogram heat map -- absolute number of voxels
        ff_histogramHeat(BarData1, BarData2, maxValue, maxValue, 50,cmapValuesHist,fov,roiName);
        numVoxels = length(BarData1);
        % axes and title
        switch fieldName
            case {'ecc'}
                if jj==1 || jj==4
                    ylabel(['pRF eccentricity for ' rm2Descript ' (deg)'])
                end
                if jj>3
                    xlabel(['pRF eccentricity for ' rm1Descript ' (deg)'])
                end
            case {'co'}
                if jj==1 % || jj==4
                    ylabel(['Variance explained for ' rm2Descript ' (%)'])
                end
                % if jj>3
                    xlabel(['Variance explained for ' rm1Descript ' (%)'])
                % end
            otherwise
                if jj==1 || jj==4
                    ylabel(['' rm2Descript ''])
                end
                if jj>3
                    xlabel(['' rm1Descript ''])
                end
        end

    end % loop over rois

    % titleName = {
    %     [strrep(ff_stringRemove(roiName, 'WangAtlas_'),'_','\_') '.' fieldName];
    % ['slope: ' num2str(meanSlope)];
    % ['ci: ' num2str(ci')];
    % [num2str(numVoxels) ' voxels']
    %      };
    % title(titleName, 'FontWeight', 'Bold');
    saveas(gcf, fullfile(sdRP,'local','png',[titleName{1} '_' fieldName '_band-2x' num2str(fov) '.png']), 'png')
    saveas(gcf, fullfile(sdRP,'local','svg',[titleName{1} '_' fieldName '_band-2x' num2str(fov) '.svg']), 'svg')
end % loop over fields

% percent above identityLine, pooled over subjects ... print out the
% percentAbovePoooled

% Percent above identity line. Bootstrapped across subjects (mixed effects)
% {   ['Percent of voxels above identity line. ']
%     ['Bootstrapped across subjects']
%     [rm2Descript ' vs. ' rm1Descript]
% }
% T = cell2table(A, 'VariableNames', {'roiName', 'fieldName', 'ciLow', 'ciHigh', 'MeanPercent'});

%% Notes
% + this is lVOTRC, do the same for V1-4,hvo1 (the ones in the paper)
% - two plots, or separate long versus short lines
% + go to white background
% - obtain numbers that show that eccc>eccw is basically noise (it will be
%      only a caption in the figure)
% + DO NOT plot any ecc diff of +- 0.5 deg
% - In the scatterplot with the light blue cones:
%    + Remove cone
%    + add +-1.5deg band
%    + below, we will only have about 3% of the data
% - compare the - and + differences, they shuold be the same for V1 v2 and
% then start changing
% - when re-running the fits, use two different HRFs (simulate that they
% will actually have size differences) and show that the effect and the
% centers will not vary
% - test: select and HRF that gives the correct size in V3 at 5deg eccen,
% and run all the analyses with this HRF
% - Use all V1, not only ventral


% test of radiality from 5 to 7 deg and 8 to 12 degs, move CB to the horizontaal and maintain the angle for words






=======
%}
