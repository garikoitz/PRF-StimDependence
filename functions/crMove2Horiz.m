function [en,an,dn] = crMove2Horiz(e,a,d)
%crMove2Horiz Moves polar coordinates to horizontal
% 
% Inputs
%   e: eccentricity
%   a: angle
%   d: by how many degrees
% 
% Outputs
%   en: eccentricity. It will be the same as e
%   an: angle
% 
% TEST: 
% Generate random e,a,d and see if we can plot the whole circle of
% posibilities
% {
%% Simulate random values in right hemifield

% Eccen
rng(12345);ew = 15*rand(10,1);
rng(23451);ec = 15*rand(10,1);
% Angle in radians
rng(34512);aw = pi*rand(10,1)-(pi/2);
rng(45123);ac = pi*rand(10,1)-(pi/2);





xx = mrvNewGraphWin('LineRadiality and Scatterplot','wide');
set(xx,'Position',[0.7 0.062 .3 .5]);

ff_polarPlot(cr.defaults.covfig.vfc); 
hold on; 
maxang = [];
minang = [];
for pp=1:length(ew)
    % Convert original to cartesian and plot in red
    % Convert to cartesian
    
    [xw, yw] = pol2cart(aw(pp), ew(pp));
    [xc, yc] = pol2cart(ac(pp), ec(pp));
    line([xw xc], [yw, yc],'Color', 'r', 'LineStyle', '-', 'LineWidth',1);
    plot(xc, yc,'or');
    
    % We know the polar angle and the ecc, rotate & translate
    % Rotate angle
    acr = 0;
    awr = aw(pp) - ac(pp);
    % Convert to cartesian and plot in blue before translation
    [xwr, ywr] = pol2cart(awr, ew(pp));
    [xcr, ycr] = pol2cart(acr, ec(pp));
    line([xwr xcr], [ywr, ycr],'Color', 'b', 'LineStyle', '-', 'LineWidth',1);
    plot(xcr, ycr,'ob');
    
    % Do the translation only in the x axis and plot in green
    xcrt  = 15;
    xcdif = xcrt - xcr;
    % if ew(pp) < ec(pp)
        
        xwrt  = xwr + xcdif;
        if xwrt > 15
            disp(pp)
        end
        
        % Calculate angles only when ecc(w) < ecc(c)
        angle = atand(ywr/(xcrt-xwrt));
        if angle > 0
            maxang = [maxang angle];
        end
        if angle < 0
            minang = [minang angle];
        end
    % else
        
    %      xwrt  = xwr + xcdif;
    %      if xwrt < xcrt
    %          xwrt = xcrt + (xcrt-xwrt);
    %      end
        % angle = -atand(ywr/(xwrt-xcrt));
    % end
    line([xwrt xcrt], [ywr, ycr],'Color', 'g', 'LineStyle', '-', 'LineWidth',2);
    plot(xcrt, ycr,'og');

    
        
    
    pause
end
                    


%}







end

