% The new defaults will not take effect if there are any open figures. To
% use them, we close all figures, and then repeat the first example.
close all;

% Defaults for activnet figures
width = 7.5;     % Width in inches
height = 8.75;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 10;      % Fontsize
lw = 1.5;      % LineWidth
msz = 16;       % MarkerSize

ac=[0.25 0.25 0.25];

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'DefaultAxesFontSize',fsz)
set(0, 'DefaultAxesFontName', 'Times')
set(0, 'DefaultFigureRenderer', 'painters',...
'DefaultAxesXcolor', ac, ...
'DefaultAxesYcolor', ac, ...
'DefaultAxesZcolor', ac,...
'DefaultTextColor', ac)
% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
