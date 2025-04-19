function vfc = ff_vfcDefault(site)

    vfc.prf_size        = true; 
    vfc.method          = 'max';         
    vfc.newfig          = true;                      
    vfc.nboot           = 50;                          
    vfc.normalizeRange  = true;              
    vfc.smoothSigma     = false;                
    vfc.cothresh        = 0.2;         
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
    vfc.gridLineWidth   = 2;
    vfc.cmapRange        = [0 pi]; % the range over which color bar
    vfc.cmapValues       = flipud(jetCmap(0,128)); % the colorbar values
    vfc.alphaValue       = ''; % 0.5;
    vfc.alphaValueDot    = ''; % 0.8;
    vfc.lineWidth        = 1.5;    

    switch site
        case 'CNI'
            vfc.fieldRange      = 15;
            vfc.eccthresh       = [.2 15];
            vfc.sigmaEffthresh  = [.2 15];   % sigma effect (sigmaMajor/sqrt(exponent))
            vfc.sigmaMajthresh  = [0 15];   % sigma major (before the exponent)
        case 'ISRAEL'
            vfc.fieldRange      = 7;
            vfc.eccthresh       = [0.2 7];
            vfc.sigmaEffthresh  = [0.2 7];   % sigma effect (sigmaMajor/sqrt(exponent))
            vfc.sigmaMajthresh  = [0 14];  % sigma major 
        case 'BCBL'
            vfc.fieldRange      = 9;
            vfc.eccthresh       = [0.2 9];
            vfc.sigmaEffthresh  = [0.2 9];   % sigma effect (sigmaMajor/sqrt(exponent))
            vfc.sigmaMajthresh  = [0 9];  % sigma major 
        otherwise
            error('Only BCBL, SNI or ISRAEL are allowed')
    end

end