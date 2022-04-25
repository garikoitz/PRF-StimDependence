function [RFcov, figHandle, all_models, weight, data] = rmPlotCoveragefromROImatfile(rm,vfc)
% [RFcov, figHandle, all_models, weight, data] = rmPlotCoveragefromROImatfile(rm,vfc)
% rmPlotCoverage - calulate the visual field coverage within an ROI
% adapted from serge's code.  I have already gotten all the rm model
% parameters and save them to an roi across subjects.  the idea here is to either load them
% or pass the rm struct to this function and then make the coverage plots
% its probably easiest to pass one rm at a time (i.e. one for each subject)
% the variable we are passing looks like this

% rm =
%      coords: [3x1633 single]
%     indices: [1633x1 double]
%        name: 'lV1_all_nw'
%     curScan: 1
%          vt: 'Gray'
%          co: [1x1633 double]
%      sigma1: [1x1633 double]
%      sigma2: [1x1633 double]
%       theta: [1x1633 double]
%        beta: [1633x3 double]
%          x0: [1x1633 double]
%          y0: [1x1633 double]
%          ph: [1x1633 double]
%         ecc: [1x1633 double]
%     session: '42111_MN'



%
% Before you run this script, you have to load 'variance explained', 'eccentricity',
% 'polar-angle' and 'prf size' into 'co', 'map', 'ph' and 'amp' fields, respectively
%
% OUTPUT
%  RFcov
% INPUT
%  prf_size:        0 = plot pRF center; 1 = use pRF size
%  fieldRange:      maximum eccentricity to plot (deg)
%  method:          'sum','max', 'clipped average', 'signed profile'
%  newfig:          make a new figure (1) or not (0). (-1 indicates don't plot
%                       anything, just return the coverage map.)
%  nboot:           the number of bootstrapping (0 means no bootstrapping)
%  normalizeRange:  if true, scale z axis to [0 1]
%  smoothSigma:     median smoothing default: 2 nearest values
%  cothresh:        threshold by variance explained in model
%  eccthresh:       2-vector ecc limits (default = [0 1.5*fieldRange])
%  nsamples:        num samples in square grid (default = 128)
%  weight:          any of {'fixed', 'parameter map', 'variance explained'} (default = fixed)
%  weightBeta:      use beta values from rmModel to weight pRFs (default = false)
%  threshByCoh:     if true, threshold by values in coherence map, not variance explained in model
%                       (note: these are often the same, but don't have to be)  %   addcenters
%  addcenters:      1 = superimpose dots for pRF centers; 0 = do not show centers
%  dualVEthresh:    use both pRF model VE and GLM fit VE as threshold
%

% 08/02 KA wrote
% 08/04 KA added bootstrapping
% 08/04 SD various mods
% 09/02 SD large rearrangments
% 09/08 JW allow superposition of pRF centers; various minor debugs
% 02/10 MB added method 'betasum'



% what's this do?
% hard to tell from the code
compVolume = false;


%% load different pRF parameters
% we have all this info already in the rm struct
% % Get coordinates for current ROI
% roi.coords   = viewGet(vw, 'roiCoords');
% roi.indices  = viewGet(vw, 'roiIndices');
% roi.name     = viewGet(vw, 'roiName');
% curScan      = viewGet(vw, 'curScan');
%
% % Get co and ph (vectors) for the current scan, within the
% % current ROI.
% vt      = vw.viewType;
% co      = rmCoordsGet(vt, rmModel,'varexp',     roi.indices);
% sigma1  = rmCoordsGet(vt, rmModel,'sigmamajor', roi.indices);
% sigma2  = rmCoordsGet(vt, rmModel,'sigmaminor', roi.indices);
% theta   = rmCoordsGet(vt, rmModel,'sigmatheta', roi.indices);
% beta    = rmCoordsGet(vt, rmModel,'beta',       roi.indices);
% x0      = rmCoordsGet(vt, rmModel,'x0',         roi.indices);
% y0      = rmCoordsGet(vt, rmModel,'y0',         roi.indices);
% clear rmModel
% already have ph and ecc
% [ph ecc] = cart2pol(x0, y0);
% as a hack

rm.co = 1- rm.rss./rm.rawrss;


co      = rm.co;
sigma1  = rm.sigma1; %.major;
sigma2  = rm.sigma2; % .minor;
sigma   = sigma1; 
theta   = rm.theta; %sigma.theta;
beta    = rm.beta;
x0      = rm.x0;
y0      = rm.y0;
ph      = rm.ph;
ecc     = rm.ecc;
[ph2, ecc2] = cart2pol(x0, y0);


% insert a flag asking if we want to y flip?
% every model has a y flip - can see when plotting coverage for a quarter
% field representation, and also when trying to define pixel location
% flip it here
% y0 = -y0;


% many things we might want to threshold by.  but not from the view
% If 'threshByCoh' is set (i.e., true), thresholding will be set based on the
% view struct's coherence field, instead of from the variance explained in
% the model.
% if vfc.threshByCoh, co = vw.co{curScan}(roi.indices); end
% if ~any(co), co = []; end
% in some sense might want to threshold the voxels before they get to this
% function


% Remove NaNs from subCo and subAmp that may be there if ROI
% includes volume voxels where there is no data.

NaNs = sum(isnan(co));
if NaNs
    fprintf('[%s]:WARNING:ROI includes voxels that have no data. These voxels are being ignored.',mfilename);
    notNaNs = ~isnan(co);
    co      = co(notNaNs);
    ph      = ph(notNaNs);
    ecc     = ecc(notNaNs);
end

% Find voxels which satisfy cothresh and eccthresh.
coIndices = co>vfc.cothresh & ...
    ecc>=vfc.eccthresh(1) & ecc<=vfc.eccthresh(2);
if ~any(coIndices)
    fprintf(1,'[%s]:No values above threshold.\n',mfilename);
    RFcov = zeros(vfc.nSamples);
    figHandle    = [];
    all_models   = [];
    weight       = [];
    data         = [];
    return
end

% also select by the mean map, if that's selected
if vfc.meanThresh > 0
    meanMapFile = fullfile(dataDir(vw), 'meanMap.mat');
    if ~exist(meanMapFile, 'file')
        warning(['A mean map threshold is specified, but no mean map ' ...
            'exists. This threshold will be ignored for now.'])
    elseif isempty(map{vw.curScan})
        warning(['A mean map threshold is specified, but no mean map ' ...
            'is computed for the current scan. ' ...
            'This threshold will be ignored for now.'])
    else
        load(meanMapFile, 'map');
        meanVals = map{vw.curScan};
        if NaNs
            meanVals = meanVals(notNaNs);
        end
        coIndices = coIndices & (meanVals > vfc.meanThresh);
    end
end

% check
if vfc.verbose
    fprintf(1,'[%s]:co-thresh:%.2f.\n',mfilename,vfc.cothresh);
    fprintf(1,'[%s]:ecc-thresh:[%.2f %.2f].\n',mfilename,vfc.eccthresh(1),vfc.eccthresh(2));
    fprintf(1,'[%s]:Number of voxels above thresh in ROI: %d (total=%d).\n',...
        mfilename,sum(coIndices),numel(coIndices));
end

% Pull out co and ph for desired pixels
subCo    = co(coIndices);
subPh    = ph(coIndices);
subEcc   = ecc(coIndices);
subSize1 = single(sigma1(coIndices));   % sigma major
subSize2 = single(sigma2(coIndices));   % sigma minor
subSize  = single(sigma(coIndices));    % sigma effective
subTheta = single(theta(coIndices));
subx0    = x0(coIndices);
suby0    = y0(coIndices);


% smooth sigma
if vfc.smoothSigma
    % vfc.smoothSigma should just really be the smoothing kernel
    if vfc.smoothSigma == 1
        
        vfc.smoothSigma = 3; % default
    end
    
    % what if we have less voxels  than smoothing kernel size?
    %* rl: add an if statement to just let n = 1 then
    if vfc.smoothSigma > length(rm.co)
        n=length(rm.co)
    else
        n = vfc.smoothSigma;
    end
    %* rl
    if length(subSize1) < vfc.smoothSigma; 
        n = 1; 
    end
    
    % check sigma1==sigma2
    if subSize1 == subSize2
        %         for every voxel
        for ii = 1:length(subSize1)
            %spatial smoothing, hard cooded to 3, not sure the units (gray
            %voxels)
            %compute nearest coords in stimulus space?
            %             get the squared distance of each voxel from all the others in
            %             stimulus space (i.e. distance between prf centers)
            %  sum of square roots is the deviation
            dev = sqrt(abs(subx0(ii) - subx0).^2 + abs(suby0(ii) - suby0).^2);
            %             sort those to closest
            [dev, ix] = sort(dev); %#ok<*ASGLU>
            %             get those up to n which is your smoothing kernel
            % and smooth those sigmas  this will have the effect of cleaning up the
            % eccentricity data which would be more useful elsewhere
            subSize1(ii) = median(subSize1(ix(1:n)));
        end
        subSize2 = subSize1;
        %         if you allowed for oval fits
    else
        for ii = 1:length(subSize1)
            %compute nearest coords
            dev = sqrt(abs(subx0(ii) - subx0).^2 + abs(suby0(ii) - suby0).^2);
            [dev, ix] = sort(dev);
            subSize1(ii) = median(subSize1(ix(1:n)));
            subSize2(ii) = median(subSize2(ix(1:n)));
        end
    end
end

% polar plot
% so the x and y positions are rederived at this point.   I'm not sure why
% not just use x and y?
% subX = single(subEcc .* cos(subPh));
% subY = single(subEcc .* sin(subPh));
% suppose we use x and y
subX=subx0;
subY=suby0;
% this fixes the failure to flip the gaussians when I flip x




% visual field indices
x = single( linspace(-vfc.fieldRange, vfc.fieldRange, vfc.nSamples) );
[X,Y] = meshgrid(x,x);

% gather this data to make accessible in the plot if you are going to use
% 'UserData'  note that here subx0 and suby0 are used not subX and subY

data.figHandle = gcf;
data.co        = co;
data.ph        = ph;
data.subCo     = subCo;
data.subPh     = subPh;
data.subEcc    = subEcc;
data.subx0     = subx0;
data.suby0     = suby0;
data.subSize1  = subSize1;  % sigma major
data.subSize2  = subSize2;  % sigma minor
data.X         = X;
data.Y         = Y;
data.subSize   = subSize;   % effective sigma
    

% For the pRF center plot, use a small constant pRF size
% not sure what happens when you set the size yourself
if vfc.prf_size==0
    %     changes all the sizes to be the same
    subSize1 = ones(size(subSize1)) * 0.1;
    subSize2 = ones(size(subSize2)) * 0.1;
    %    replaces ovals with circles if there were any
    subTheta = zeros(size(subTheta));
end


% various ways to weight the data when taking averages
switch lower(vfc.weight)
    case 'fixed'
        weight = ones(size(subCo));
        
    case 'parameter map'
        weight = getCurDataROI(vw,'map',curScan,roi.coords);
        weight = weight(coIndices);
        
    case {'variance explained', 'varexp', 've'}
        weight = subCo;
        
    otherwise
        error('Unknown weight parameter: %s',vfc.weight);
end


if vfc.weightBeta==1
    weight = weight .* beta(coIndices);
end
weight = single(weight);

%% special case: for the 'density' coverage option, we don't need to
%% do a lot of memory-hungry steps like making all pRFs. So, I've set those
%% computations aside in their own subroutine. (ras)
if isequal( lower(vfc.method), 'density' )
    RFcov = prfCoverageDensityMap(vw, subx0, suby0, subSize1, X, Y);
    
    all_models = []; % not created for this option
    if vfc.newfig==-1
        figHandle = [];
    else
        figHandle = createCoveragePlot(vw, RFcov, vfc, roi, data, m);
    end
    
    return
end


%% make all pRFs:
% make in small steps so we don't go into swap space for large ROIs
% number of pRFs
n = numel(subX);
% s is step size for pRF calculation.  since it makes many prfs
% simultaneously you have to put an upper limit on that which is apparently
% 1000
s = [(1:ceil(n./1000):n-2) n+1];

% NOTE: (RL)
% For the line above (which assumes that we have at least 3 voxels,
% probably for median smoothing),  s is an empty vector when n is less than 
% 3 voxels -- and an empty rfcov is  returned when we try to plot the 
% coverage. So modify s accordingly for these edge cases: 
if n < 3
    s = [1:n+1]; 
end

% matrix for storing models image size (as vector) x num of prfs
all_models = zeros( numel(X), n, 'single' );
fprintf(1,'[%s]:Making %d pRFs:...', mfilename, n);
drawnow;
% in here is the where the rfs are made.  in case we want them to swap
% sides
% so for each step up to the last one
for n=1:numel(s)-1,
    % make rfs
    rf   = rfGaussian2d(X(:), Y(:),...
        subSize1(s(n):s(n+1)-1), ...
        subSize2(s(n):s(n+1)-1), ...
        subTheta(s(n):s(n+1)-1), ...
        subX(s(n):s(n+1)-1), ...
        subY(s(n):s(n+1)-1));
    all_models(:,s(n):s(n+1)-1) = rf;
end;

n; % debugging purposes, remove anytime 

% so rf is a 2d matrix
% each column is a centered prf fit
% as many rows as there will be pixels in the image (so X^2)
% these columns are stacked together in all_models


clear n s rf pred;
fprintf(1, 'Done.\n');
drawnow;

% Correct volume
% not sure what this is doing for us
if compVolume
    tmp = ones(size(all_models, 1), 1, 'single');
    %     square the size of all your sigmas
    vol = sigma1(coIndices).^2;
    %     then multiply them by 2 pi
    vol = vol * (2 * pi);
    
    all_models = all_models ./ (tmp * vol);
end

% For the pRF center plot, put a constant value (1) within each Gaussian
if vfc.prf_size==0
    all_models(all_models>0.1)=1;
end

% weight all models
if isequal( lower(vfc.weight), 'fixed' )
    % if the weights are even, we avoid the redundant, memory-hungry
    % multiplication step that would otherwise be done.
    all_models_weighted = all_models;
else
    tmp = ones(size(all_models, 1), 1, 'single');
    all_models_weighted = all_models .* (tmp * weight);
    clear tmp
end

%% Different ways of combining them:
% 1) bootstrap (yes, no) 2) which statistic (sum, max, etc),
% bootstrap

% If we are only working with 1 voxel, bootstrapping will not make a
% difference. We turn it off because the bootstp function does not handle
% this case well. 
if size(subX,2) == 1
    vfc.nboot = 0; 
end

if vfc.nboot>0
    if isempty(which('bootstrp'))
        warndlg('Bootstrap requires statistics toolbox');
        RFcov = [];
        return;
    end
    all_models(isnan(all_models))=0;
    
    switch lower(vfc.method)
        %         sample with replacement. compute average of all pRFs
        case {'sum','add','avg','average everything'}
            m = bootstrp(vfc.nboot, @mean, all_models');
            %         sample with replacement.  take max value at each location
        case {'max','profile','maximum profile' 'maximum'}
            m = bootstrp(vfc.nboot, @max, all_models');
            
        otherwise
            error('Unknown method %s',vfc.method)
    end
    % take mean of bootstraps
    RFcov=mean(m,1)';
    
    % no bootstrap
else
    
    m = []; 
    switch lower(vfc.method)
        
        % coverage = sum(pRF(i)*w(i)) / (sum(pRF(i))
        case {'beta-sum','betasum','weight average'}
            RFcov = sum(all_models_weighted, 2) ./ sum(all_models,2);
            
            % coverage = sum(pRF(i)*w(i)) / (sum(pRF(i)) + clipping
        case {'clipped beta-sum','clippedbeta','clipped weight average'}
            % set all pRF beyond 2 sigmas to zero
            clipval = exp( -.5 *((2./1).^2));
            all_models(all_models<clipval) = 0;
            n = all_models > 0;
            
            % recompute all_models_weighted
            tmp = ones( size(all_models,1), 1, 'single' );
            all_models_weighted = all_models .* (tmp*weight);
            
            % compute weighted clipped sum/average
            sumn = sum(n,2);
            mask = sumn==0;
            sumn(mask) = 1; % prevent dividing by 0
            RFcov = sum(all_models_weighted,2) ./ sum(all_models,2);
            RFcov(mask) = 0;
            
            %clip to zero if n<clipn
            if isnumeric(vfc.clipn)
                RFcov(sumn<=vfc.clipn) = 0;
            end
            
            % coverage = sum(pRF(i)*w(i)) / (sum(w(i))
        case {'sum','add','avg','average','prf average'}
            RFcov = sum(all_models_weighted, 2) ./ sum(weight);
            
            % coverage = sum(pRF(i)*w(i)) / (sum(w(i)) + clipping
        case {'clipped average','clipped','clipped prf average'}
            % set all pRF beyond 2 sigmas to zero
            clipval = exp( -.5 *((2./1).^2));
            all_models(all_models<clipval) = 0;
            n = all_models > 0;
            
            % recompute all_models_weighted
            tmp = ones( size(all_models,1), 1, 'single' );
            all_models_weighted = all_models .* (tmp*weight);
            
            % compute weighted clipped mean
            sumn = sum(weight.*n);
            mask = sumn==0;
            sumn(mask) = 1; % prevent dividing by 0
            RFcov = sum(all_models_weighted,2) ./ sumn;
            RFcov(mask) = 0;
            
            %clip to zero if n<clipn
            if isnumeric(vfc.clipn)
                RFcov(sumn<=vfc.clipn) = 0;
            end
            
            % coverage = max(pRF(i))
        case {'maximum profile', 'max', 'maximum'}
            RFcov = max(all_models_weighted,[],2);
            
        case {'signed profile'}
            RFcov  = max(all_models_weighted,[],2);
            covmin = min(all_models_weighted,[],2);
            ii = RFcov<abs(covmin);
            RFcov(ii)=covmin(ii);
            
        case {'p','probability','weighted statistic corrected for upsampling'}
            RFcov = zeros(vfc.nSamples);
            
            % I guess this upsample factor assumes your functional data are
            % 2.5 x 2.5 x 3 mm?
            upsamplefactor = 2.5*2.5*3; % sigh.....
            for ii = 1:size(all_models,1)
                s = wstat(all_models(ii,:),weight,upsamplefactor);
                if isfinite(s.tval)
                    RFcov(ii) = 1 - t2p(s.tval,1,s.df);
                end
            end
            
        otherwise
            error('Unknown method %s',vfc.method)
    end
end

% convert 1D to 2D
RFcov = reshape( RFcov, [1 1] .* sqrt(numel(RFcov)) );

% When no voxels exceed threshold, return nan matrix rather than empty
% matrix
if sum(size(RFcov))==0
    RFcov=nan(nSamples,nSamples);
end

% if the newfig flag is set to -1, just return the image
if vfc.newfig==-1,
    figHandle = [];
else
    figHandle = createCoveragePlot( RFcov, vfc,rm, data, m);
end


return
% /--------------------------------------------------------------------/ %




% /--------------------------------------------------------------------/ %
function figHandle = createCoveragePlot(RFcov, vfc,rm, data, m)
% plotting subroutine for rmPlotCoverage. Broken off by ras 10/2009.
% This is not working for me, can't write figures in loop
% if vfc.newfig
%     figHandle = figure('Color', 'w');
% else
%     figHandle = selectGraphWin;
% end

figHandle = mrvNewGraphWin;


% headerStr = sprintf('Visual field coverage, ROI %s, scan %i', ...
% 					roi.name, vw.curScan);
set(gcf, 'Name', rm.session);


% normalize the color plots to 1
if vfc.normalizeRange,
    rfMax = max(RFcov(:));
else
    rfMax = 1;
end

img = RFcov ./ rfMax;
mask = makecircle(length(img));
img = img .* mask;
imagesc(data.X(1,:), data.Y(:,1), img);

set(gca, 'YDir', 'normal');
grid on

colormap(vfc.cmap);
colorbar;

% start plotting
hold on;

% t = 0:.01:2*pi;
%
% % rings every 5 deg
% for n=(1:3)/3*vfc.fieldRange
%     polar(t,ones(size(t))*n,'w');
% end
% plot([0 0],[-vfc.fieldRange vfc.fieldRange],'w')
% plot([-sqrt(vfc.fieldRange^2/2) sqrt(vfc.fieldRange^2/2)],[-sqrt(vfc.fieldRange^2/2) sqrt(vfc.fieldRange^2/2)],'w')
% plot([-vfc.fieldRange vfc.fieldRange],[0 0],'w')
% plot([-sqrt(vfc.fieldRange^2/2) sqrt(vfc.fieldRange^2/2)],[sqrt(vfc.fieldRange^2/2) -sqrt(vfc.fieldRange^2/2)],'w')


% add polar grid on top
p.ringTicks = (1:3)/3*vfc.fieldRange;
p.color = 'w';
if isfield(vfc, 'backgroundColor')
    p.backgroundColor = vfc.backgroundColor;
end
if isfield(vfc, 'tickLabel')
    p.tickLabel = vfc.tickLabel;
end

polarPlot([], p);

% add pRF centers if requested
if vfc.addCenters,
    inds = data.subEcc < vfc.fieldRange;
    plot(data.subx0(inds), data.suby0(inds), '.', ...
        'Color', [.5 .5 .5], 'MarkerSize', 4);
end

if isfield(vfc, 'ellipsePlot') || isfield(vfc, 'contourPlot')
    % get the x and y points of the specified contour level
    [contourMatrix, contourCoordsX, contourCoordsY] = ...
        ff_contourMatrix_makeFromMatrix(RFcov,vfc,vfc.contourLevel);

    % transform so that we can plot it on the polar plot
    contourX = contourCoordsX/vfc.nSamples*(2*vfc.fieldRange) - vfc.fieldRange; 
    contourY = contourCoordsY/vfc.nSamples*(2*vfc.fieldRange) - vfc.fieldRange;  
end

if isfield(vfc, 'ellipsePlot')
    if vfc.ellipsePlot
        ellipse_t = ff_fit_ellipse(contourX, contourY, gca);
    end
end

if isfield(vfc, 'contourPlot')
    if vfc.contourPlot
        plot(contourX, contourY, '--', 'Color', vfc.contourColor, 'MarkerSize', 2, 'LineWidth',1.5);
    end
end

if isfield(vfc, 'contourBootstrap')
    if vfc.contourBootstrap
        % draws  ellipses or contours bootstrapped vfc.nboot times depending on
        % whether whether each of them are true. assumes it will draw it on the
        % current figure
        ff_contourBootstrapDraw(RFcov,vfc, data, m);
    end
end



% scale z-axis
if vfc.normalizeRange
    if isequal( lower(vfc.method), 'maximum profile' )
        caxis([.5 1]);
    else
        caxis([0 1]);
    end
else
    if min(RFcov(:))>=0
        caxis([0 ceil(max(RFcov(:)))]);
    else
        caxis([-1 1] * ceil(max(abs(RFcov(:)))));
    end
end
axis image;   % axis square;
xlim([-vfc.fieldRange vfc.fieldRange])
ylim([-vfc.fieldRange vfc.fieldRange])

title(rm.name, 'FontSize', 24, 'Interpreter', 'none');


% Save the data in gca('UserData')
set(gca, 'UserData', data);

return;
% /------------------------------------------------------------------/ %




% /------------------------------------------------------------------/ %
function RFcov = prfCoverageDensityMap(vw, x0, y0, sigma, X, Y) %#ok<INUSL>
% for each point (x, y) in visual space, this returns
% the proportion of voxels in the ROI for which (x, y) is
% within one standard deviation of the pRF center.
mask = NaN( size(X, 1), size(X, 2), length(x0) );

for v = 1:length(x0)
    % make a binary mask within one sigma of the center
    R = sqrt( (X - x0(v)) .^ 2 + (Y - y0(v)) .^ 2 );
    mask(:,:,v) = ( R < 2*sigma(v) );
end

% average (sum?) across all masks
RFcov = nansum(mask, 3);

return
% /------------------------------------------------------------------/ %




% /------------------------------------------------------------------/ %
function vfc = rmPlotCoverageDialog(vfc)
%% dialog to get parameters for rmPlotCoverage.
dlg(1).fieldName = 'method';
dlg(end).style = 'popup';
dlg(end).list = {'maximum profile' 'average' 'clipped average' ...
    'density' 'signed profile' 'probability'};
dlg(end).string = 'Method for combining pRFs?';
dlg(end).value = vfc.(dlg(end).fieldName);


dlg(end+1).fieldName = 'weight';
dlg(end).style = 'popup';
dlg(end).list = {'fixed' 'variance explained' 'parameter map'};
dlg(end).string = 'Method for weighting pRFs?';
dlg(end).value = vfc.(dlg(end).fieldName);


dlg(end+1).fieldName = 'fieldRange';
dlg(end).style = 'number';
dlg(end).string = 'Visual Field Range (deg)?';
dlg(end).value = vfc.(dlg(end).fieldName);


dlg(end+1).fieldName = 'nboot';
dlg(end).style = 'number';
dlg(end).string = 'Number of bootstrapping steps?';
dlg(end).value = vfc.(dlg(end).fieldName);


dlg(end+1).fieldName = 'cmap';
dlg(end).style = 'popup';
dlg(end).list = mrvColorMaps;
dlg(end).string = 'If plotting, color map for coverage?';
dlg(end).value = vfc.(dlg(end).fieldName);


dlg(end+1).fieldName = 'normalizeRange';
dlg(end).style = 'checkbox';
dlg(end).list = {};
dlg(end).string = 'Normalize data range to [0 1]';
dlg(end).value = vfc.(dlg(end).fieldName);


dlg(end+1).fieldName = 'smoothSigma';
dlg(end).style = 'checkbox';
dlg(end).list = {};
dlg(end).string = 'Smooth sigma (medianfilter)';
dlg(end).value = vfc.(dlg(end).fieldName);


dlg(end+1).fieldName = 'prf_size';
dlg(end).style = 'checkbox';
dlg(end).list = {};
dlg(end).string = 'Use pRF sizes from model';
dlg(end).value = vfc.(dlg(end).fieldName);


dlg(end+1).fieldName = 'newfig';
dlg(end).style = 'checkbox';
dlg(end).string = 'Show results in new figure';
dlg(end).value = vfc.(dlg(end).fieldName);


dlg(end+1).fieldName = 'addCenters';
dlg(end).style = 'checkbox';
dlg(end).string = 'Add dots to show pRF centers';
dlg(end).value = vfc.(dlg(end).fieldName);


[resp ok] = generalDialog(dlg, mfilename);
if ~ok
    error('User Aborted.')
end
drawnow;

vfc = mergeStructures(vfc, resp);

return
