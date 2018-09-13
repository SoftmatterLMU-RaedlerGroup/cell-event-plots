%plotCorrelation_single plots correlation scatter plots.
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
	save_path = fullfile(out_dir, [ getTime '_' out_token_temp 'eventTimeCorr_singleMeasurements.ps' ]);
else
	save_path = false;
end

%% Plot combinations in separate subplots
fhc = cell(1, height(corrT));

for i_cmb = 1:height(corrT)
	%% Get measurements of files/traces
	file_map = corrT.file_map{i_cmb,1};
	files = zeros(size(file_map, 1), 1);

	for i_fl = 1:size(file_map, 1)
		files(file_map(i_fl, 2):file_map(i_fl, 3)) = file_map(i_fl, 1);
	end

	[~, files_idcs] = ismember(files, fileT.index_f);
	[measurements, ~, measurements_idcs] = ...
		unique(fileT.measurement(files_idcs));
	measurement_colors = lines(length(measurements));
	measurement_symbols = cell(size(measurements, 1), 1);

	%% Perform plot
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
			corrT.markers{i_cmb, 2}), ...
			'Interpreter', 'tex' ...
		)

	% Plot the event times
	if isempty(corrT.t_event{i_cmb})
		% Skip empty combinations (to avoid matrix index error)
		warning('No data found for combination %s/%s.', ...
			corrT.markers{i_cmb, 1}, ...
			corrT.markers{i_cmb, 2});
		%continue

	else
		% Create a scatter plot
		for i_trace = 1:length(files)
			this_data = zeros(1,2);
			doNotPlot = true;

			if any(corrT.found{i_cmb}(i_trace,[1,2]))
				this_data(1) = corrT.t_event{i_cmb}(i_trace, 1);
				doNotPlot = false;
			end
			if any(corrT.found{i_cmb}(i_trace,[1,3]))
				this_data(2) = corrT.t_event{i_cmb}(i_trace, 2);
				doNotPlot = false;
			end

			% Do not plot traces with neither event found
			if doNotPlot
				continue
			end

			% Find symbol for plotting
			% "." for both events found, "+" for only one event found
			if corrT.found{i_cmb}(i_trace,1)
				plt_symb = '.';
			else
				plt_symb = '+';
			end

			% Get measurement index
			i_meas = measurements_idcs(i_trace);

			% Plot
			p = plot(ax, this_data(1), this_data(2), ...
				plt_symb, 'Color', measurement_colors(i_meas,:));

			% Add symbol to legend
			ms = measurement_symbols{i_meas};
			if isempty(ms) || ( (ms.Marker == '+') && (plt_symb == '.') )
				measurement_symbols{i_meas} = p;
			end
		end
	end

	lgd = legend(ax, [measurement_symbols{:}], ...
		{measurements{:}}, 'Location', 'eastoutside');
	title(lgd, 'Measurements')

	%% Export figure if desired
	if save_path
		fhc{i_cmb}.PaperType = 'a5';
		fhc{i_cmb}.PaperOrientation = 'landscape';
		fhc{i_cmb}.PaperPositionMode = 'manual';
		fhc{i_cmb}.PaperUnits = 'normalized';
		fhc{i_cmb}.PaperPosition = [0 0 1 1];
		fhc{i_cmb}.PaperUnits = 'centimeters';

		print(fhc{i_cmb}, '-painters', save_path, '-dpsc', '-append')
	end
end

%% Cleanup
clear t_max save_path out_token_temp conditions text_arrow i_meas
clear i_cmb ax l this_data doNotPlot colors_idcs file_map files files_idcs
clear i_fl i_trace measurement_colors plt_symb lgd p ms
