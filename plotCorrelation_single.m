%plotCorrelation_single plots correlation scatter plots.
%
% Set `makeEllipses` to true to plot error ellipses
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

if ~exist('t_max', 'var')
	t_max = 30;
end

%% Get filename for output file
% Check output directory
if exist('out_dir', 'var')
	% Ensure target directory exists
	if ~isempty(out_dir) && ~exist(out_dir, 'dir')
		[~,~,~] = mkdir(out_dir);
	end

	% Check conditions
	conditions = unique(fileT.condition(ismember(uint32(fileT.index_f), ...
		[corrT.file_map{1,1}(:,1); corrT.file_map{1,2}(:,1)])));

	if ~isempty(conditions)
		out_token_temp = [ strjoin(conditions, '-') '_' ];
		if isstring(out_token_temp)
			out_token_temp = [ out_token_temp{:} ];
		end
	else
		out_token_temp = '';
	end

	% Assemble file name
	save_path = fullfile(out_dir, [ getTime '_' out_token_temp 'eventTimeCorr_single.ps' ]);
else
	save_path = false;
end

%% Check whether to create error ellipses
if exist('makeEllipses', 'var')
	if ~islogical(makeEllipses)
		makeEllipses = logical(makeEllipses);
	end
else
	makeEllipses = false;
end

if makeEllipses
	ellipseCenters = NaN(height(corrT), 2);
else
	ellipseCenters = [];
end
	

%% Plot combinations in separate subplots
fhc = cell(1, height(corrT));

for i_cmb = 1:height(corrT)
	% Prepare figure and axes
	fhc{i_cmb} = figure;
	ax = axes(fhc{i_cmb});
	hold(ax, 'on')

	% Formating and labeling axis
	ax.Box = 'On';
	ax.XLim = [0 t_max];
	ax.YLim = [0 t_max];

	l = line(ax, xlim, ylim, 'Color', 'k');
	uistack(l, 'bottom')

	xlabel(ax, ['$t_\mathrm{event}(\mathrm{' corrT.markers{i_cmb,1} '})$ [h]'], 'interpreter','latex')
	ylabel(ax, ['$t_\mathrm{event}(\mathrm{' corrT.markers{i_cmb,2} '})$ [h]'], 'interpreter','latex')

	text_arrow = char(8594);
	title(ax, ...
		sprintf( '\\color[rgb]{%f,%f,%f}%s: %s %s %s ', ...
			corrT.color(i_cmb,:), ...
			strjoin(conditions, '/'), ...
			corrT.markers{i_cmb, 1}, ...
			text_arrow, ...
			corrT.markers{i_cmb, 2}) ...
		)

	% Plot the event times
	this_data = corrT.t_event{i_cmb};

	if isempty(this_data)
		% Skip empty combinations (to avoid matrix index error)
		warning('No data found for combination %s/%s.', ...
			corrT.markers{i_cmb, 1}, ...
			corrT.markers{i_cmb, 2});
		%continue

	else
		% Create a scatter plot
		plot(ax, this_data(corrT.found{i_cmb}(:,1),1), this_data(corrT.found{i_cmb}(:,1),2), '.', 'Color', corrT.color(i_cmb,:))
		plot(ax, this_data(corrT.found{i_cmb}(:,2),1), zeros(sum(corrT.found{i_cmb}(:,2)),1), '+', 'Color', corrT.color(i_cmb,:))
		plot(ax, zeros(sum(corrT.found{i_cmb}(:,3)),1), this_data(corrT.found{i_cmb}(:,3),2), '+', 'Color', corrT.color(i_cmb,:))

		% Plot error ellipse
		if makeEllipses
			[ellips, comX, comY] = error_ellipse(this_data(corrT.found{i_cmb}(:,1),:));
			ellips = ellips + [comX comY];
			plot(ax, ellips(:,1), ellips(:,2), '-k');
			plot(ax, comX, comY, '+k');

			ellipseCenters(i_cmb,:) = [comX, comY];
		end
	end

	%% Export figure if desired
	if save_path
		set(fhc{i_cmb}, ...
			'PaperPositionMode', 'manual', ...
			'PaperUnits', 'centimeters', ...
			'PaperSize', [17 15], ...
			'PaperPosition', [0 0 17 15] ...
			)
		print(fhc{i_cmb}, '-painters', save_path, '-dpsc', '-append')
	end
end

%% Cleanup
clear t_max save_path out_token_temp conditions
clear i_cmb ax l this_data text_arrow