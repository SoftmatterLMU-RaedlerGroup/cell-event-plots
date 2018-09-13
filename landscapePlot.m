%landscapePlot plots event time distributions as 3-dimensional log-normal distributions
%
% The user is prompted which combination to plot and how many log-normal
% distributions to fit against the data.
% The figure can be saved to a file by rightclicking (outside of the axes)
% and selecting the print quality (vector or raster graphics).
%
% Copyright © 2018 Daniel Woschée <daniel.woschee@physik.lmu.de>
% Faculty of Physics / Ludwig-Maximilians-Universität München
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License version 2
% as published by the Free Software Foundation.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%%

%% Find out which combination to plot
% Prompt user for combinations to analyze
[comb_sel, isOK] = listdlg('ListString', join(corrT.markers, ' --> '), ...
	'Name', 'Combination selection', ...
	'PromptString', 'Select the combinations you want to load:', ...
	'SelectionMode', 'single');
if ~isOK || isempty(comb_sel)
	return
end

% Prompt user for number of probability distributions
[n_dists, isOK] = listdlg('ListString', string(1:5), ...
	'Name', 'Number of distributions', ...
	'PromptString', 'Please select the number of probability distributions to use:', ...
	'SelectionMode', 'single');
if ~isOK || isempty(n_dists)
	return
end

%% Assess conditions of plotted traces
% Extract conditions from `fileT`
conditions = unique(fileT.condition(ismember(fileT.index_f, ...
	[corrT.file_map{comb_sel,1}(:,1); corrT.file_map{comb_sel,2}(:,1)])));

% Create condition token for use in output file name
% (must be char vector because `fullfile` cannot cope with string input
condition_fn = strjoin(conditions, '-');
if strlength(condition_fn) > 0
	condition_fn = [condition_fn{1} '_'];
else
	condition_fn = '';
end

% Create condition token for use in figure title
condition_ttl = strjoin(conditions, ', ');
if strlength(condition_ttl) > 0
	condition_ttl = [condition_ttl{1} ': '];
end


%% Perform fit
paramBounds = [0, -1+eps, 0, 0, 0, 0; 1, +1-eps, 4, 4, 10, 10];
paramBounds = repmat(paramBounds, 1, n_dists);

t_x(1) = eps;

if ~exist('t_max', 'var')
	t_max = 30;
end

% 3-dim histogram
[hist_count, hist_centers] = hist3count( ...
	corrT.t_event{comb_sel}(corrT.found{comb_sel}(:,1),1), ...
	corrT.t_event{comb_sel}(corrT.found{comb_sel}(:,1),2), ...
	30, 30);

fit_x = hist_centers.x;
fit_y = hist_centers.y;
if ~isrow(fit_x)
	fit_x = fit_x.';
end
if ~iscolumn(fit_y)
	fit_y = fit_y.';
end

fit_coord = zeros(length(fit_y), length(fit_x), 2);
fit_coord(:, :, 1) = repmat(fit_x, length(fit_y), 1);
fit_coord(:, :, 2) = repmat(fit_y, 1, length(fit_x));

fitPar = fit_mle( ...
	@(coord,params)logNormBiN(n_dists, coord, params), ...
	fit_coord, ...
	hist_count.matrix / hist_count.normfac, ...
	'ParameterBounds', paramBounds, ...
	'ParameterList', 1:size(paramBounds, 2) ...
	);

%% Plot data
% Get time vectors for plot
plt_res = 500;
t_plt = struct();

t_plt.x = linspace(0, t_max, plt_res)'; t_plt.x(t_plt.x == 0) = eps;
t_plt.y = linspace(0, t_max, plt_res)'; t_plt.y(t_plt.y == 0) = eps;

% Write to figure
fh = figure;
ax = axes(fh);

idx_count_nonzero = hist_count.linearized ~= 0;
stem3(ax, hist_centers.linearized(idx_count_nonzero,1), ...
	hist_centers.linearized(idx_count_nonzero,2), ...
	hist_count.linearized(idx_count_nonzero) / hist_count.normfac, '.r')
% plot3(ax, hist_centers.linearized(idx_count_nonzero,1), ...
% 	hist_centers.linearized(idx_count_nonzero,2), ...
% 	hist_count.linearized(idx_count_nonzero) / hist_count.normfac, '.r')
hold on

surface(ax, t_plt.x, t_plt.y, logNormBiN(n_dists, {t_plt.x.', t_plt.y}, fitPar))
shading interp
xlabel(ax, sprintf('$t_\\mathrm{event}(\\mathrm{%s})$ [h]', ...
	corrT.markers{comb_sel,1}), 'Interpreter', 'latex')
ylabel(ax, sprintf('$t_\\mathrm{event}(\\mathrm{%s})$ [h]', ...
	corrT.markers{comb_sel,2}), 'Interpreter', 'latex')
zlabel(ax, 'Probability')
title(ax, sprintf('%sDistribution of event times for %s/%s', ...
	condition_ttl, corrT.markers{comb_sel,:}))

%% Prepare figure for plotting
% Define size of output figure
fh.Units = 'centimeters';
fh.PaperPositionMode = 'manual';
fh.PaperSize = [21 17];
fh.PaperPosition = [0 0 fh.PaperSize];

% Set information needed for exporting image to file
fh.UserData.corrT = corrT(comb_sel,:);
fh.UserData.condition_fn = condition_fn;
fh.UserData.out_dir = '';
if exist('out_dir', 'var')
	fh.UserData.out_dir = out_dir;
end

% Build context menu for plotting the figure in desired quality
cm = uicontextmenu;
uimenu(cm, 'Label', 'Print as vector graphics', ...
	'Callback', {@plotThis, 'vector'});
uimenu(cm, 'Label', 'Print as raster graphics', ...
	'Callback', {@plotThis, 'raster'});
fh.UIContextMenu = cm;

%% Callback function for plotting context menu
function plotThis(gh, ~, quality)
% Get requested quality
switch quality
	case 'vector'
		suffix = '.pdf';
		format = '-dpdf';
		renderer = "-painters";

	case 'raster'
		suffix = '.png';
		format = '-dpng';
		renderer = ["-opengl" "-r1000"];

	otherwise
		return
end

% Get figure handle
while ~isa(gh, 'matlab.ui.Figure')
	gh = gh.Parent;
end

% Assess target directory
if  ~exist(gh.UserData.out_dir, 'dir')
	out_dir = uigetdir('', 'Select directory to save figure');
	if out_dir == 0
		return
	else
		gh.UserData.out_dir = out_dir;
	end
end

% Get time as unique file identifier
time_now = clock;
time_now = sprintf('%4d-%02d-%02d–%02d%02d', ...
	time_now(1), time_now(2), time_now(3), time_now(4), time_now(5));

% Assemble file name
save_path = fullfile(gh.UserData.out_dir, [ time_now '_landscape_' ...
	gh.UserData.condition_fn char(join(gh.UserData.corrT.markers, '-')) suffix ]);

% Save the figure to the given file
print(gh, renderer{:}, save_path, format)
fprintf('Figure written to:\n%s\n', save_path);

end
