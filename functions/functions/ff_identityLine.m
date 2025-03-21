function p = ff_identityLine(ax, lineColor)
%Draw an identity line on the current axis
%
%   p = identityLine(ax, lineColor)
%
% Example:
%   plot(1:10,randn(1,10),'o')
%   identityLine(gca, [.5 .5 .5]);
%
% (c) Stanford VISTA Team

if notDefined('ax'), ax = gca; end

% Minimum and maximum of axes
xlim = get(ax,'xlim');
ylim = get(ax,'ylim');
mn = min(xlim(1),ylim(1));
mx = max(xlim(2),ylim(2));

% Here's the line
p = line([mn mx],[mn mx],'color',lineColor,'linestyle','--');

% Set line properties.  These probably want to come in as an argument
set(p,'linewidth',2);
grid on

return