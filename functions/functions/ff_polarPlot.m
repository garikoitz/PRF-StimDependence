function ff_polarPlot(vfc)
% ff_polarPlot(vfc)
%
% add polar grid on top

p.ringTicks = (1:3)/3*vfc.fieldRange;
p.gridColor = vfc.gridColor;                % default is [.6 .6 .6]
p.backgroundColor = vfc.backgroundColor;    % what it will be displayed as (fig)
p.gridLineWidth=vfc.gridLineWidth;

% what the png will be saved as
if ~isfield(vfc, 'fillColor')
    p.fillColor = [1 1 1];
else
    p.fillColor = vfc.fillColor;
end

                

if ~isfield(vfc, 'tickLabel')
    p.tickLabel = true;
else
    p.tickLabel = vfc.tickLabel;
end

polarPlot([], p);

end