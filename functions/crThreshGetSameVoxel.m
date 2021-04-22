function [R,C_data,cr]=crThreshGetSameVoxel(cr,...
                                         rmroiCell, ...
                                         list_subInds,...
                                         list_roiNames,...
                                         varargin)
%% Varargin
% Make varargin lower case, remove white spaces...
varargin = mrvParamFormat(varargin);
% Parse
p = inputParser;
p.addRequired('cr'            , @isstruct);
p.addRequired('rmroicell'     , @iscell);
p.addRequired('list_subInds'  , @isnumeric);
p.addRequired('list_roiNames' , @iscell);

p.addOptional('cothresh'      , 0.2,  @isnumeric);
p.addOptional('fieldrange'    , 15 ,  @isnumeric);

% Parse. Assign result inside each case
p.parse(cr, rmroiCell, list_subInds, list_roiNames, varargin{:});
% Read here only the generic ones
cothresh      = p.Results.cothresh;
fieldRange    = p.Results.fieldrange;

                                     
                                     
                                     
                                     
%%                                     
%  Threshold and get identical voxels for each subject
% values to threshold the RM struct by
% vfc = ff_vfcDefault_Hebrew;
% The defaults are different for the two projects, so for now use the most
% restrictive one
% vfc threshold
vfc.prf_size         = true;
vfc.fieldRange       = fieldRange; % 15 CNI % 7 Israel;
vfc.method           = 'max';
vfc.newfig           = true;
vfc.nboot            = 50;
vfc.normalizeRange   = true;
vfc.smoothSigma      = false;
% Thresholds
vfc.sigmaEffthresh   = [.2 15];  % [0  7]; % sigma effect (sigmaMajor/sqrt(exponent))
vfc.sigmaMajthresh   = [0.5 8];   % [0 14]; % sigma major (before the exponent)
vfc.cothresh         = cothresh; % Variance explained, 20%
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
L_data  = cell(1, numRois);
X_rm1   = cell(1, numRois);
Y_rm1   = cell(1, numRois);
X_rm2   = cell(1, numRois);
Y_rm2   = cell(1, numRois);
Ecc_rm1 = cell(1, numRois);
Ecc_rm2 = cell(1, numRois);

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


% Prepare the output
R = struct();
R.rmroiCellSameVox = rmroiCellSameVox;
R.L_data = L_data;
R.X_rm1 = X_rm1;
R.Y_rm1 = Y_rm1;
R.X_rm2 = X_rm2;
R.Y_rm2 = Y_rm2;
R.Ecc_rm1 = Ecc_rm1;
R.Ecc_rm2 = Ecc_rm2;
R.Ph_rm1 = Ph_rm1; 
R.Ph_rm2 = Ph_rm2;
R.Sm_rm1 = Sm_rm1; 
R.Sm_rm2 = Sm_rm2;

end
