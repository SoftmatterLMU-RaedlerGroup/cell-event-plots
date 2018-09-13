% plot_eventTimeCorrelations.m
% This script plots graphs of correlated event times for different markers.
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

%% Define Parameters
index_tevent = 2;			% index of event time in tracesA
time_scale = 1;			% time scaling factor (for conversion to hours)
% marker_combinations = cell(0,2);	% n-by-2 cellstring of marker tokens
% z.B. marker_combinations = { {'lyso','ros'}, {'psiva','pi'} };
marker_combinations = { ...
	{'tmrm','ros'}, {'lyso','tmrm'}, {'lyso','ros'}; ...
	{'psiva','tmrm'}, {'psiva','pi'}, {'lyso','pi'} };
t_max = 30;	% Maximum time to be shown

% Define white list of conditions to be used
if ~exist('restrict_to_conditions', 'var')
	restrict_to_conditions = {};
end

% Convert t_max to hours
t_max = t_max * time_scale;

% Tolerance for event time equality (in [h])
thresh_eq = 1/6;

% Define colors to be used (must have same size as marker_combinations
combination_colors = { ...
	[0.2275, 0.7765, 0.9412], ... % MOMP–ROS
	[0.2314, 0.4118, 0.6980], ... % LMP–MOMP
	[0.2000, 0.1569, 0.4000]; ... % LMP–ROS
	[0.9373, 0.2314, 0.1373], ... % PhS–MOMP
	[0.6118, 0.0824, 0.0510], ... % PhS–PMP
	[0.2275, 0.0196, 0.0196] ... % LMP–PMP
	};
combination_colors = mat2cell(lines(6), ones(1,6), 3);

% Define additional token to be added to output file names
if ~exist('out_token', 'var')
	out_token = [];
end

%% Retrieve debug status
if ~exist('DEBUG', 'var')
	DEBUG = false;
end

% Log file
if DEBUG
	fid = fopen('out/log.txt','w');
end

%% Import data
exist_files = exist('fileS', 'var') && exist('out_dir', 'var');
exist_data = exist('tracesA', 'var') && exist('metaS', 'var') ...
	&& exist('dataF', 'var');
if ~exist_data
	if ~exist_files
		[fileT, out_dir, fileS] = getFiles;
	end
	[tracesA, metaS, dataF] = readinMulticolorData(fileS);
end

%% Search relevant traces
% Allocate variables for relevant traces
ncombs = numel(marker_combinations);
data_files = cell(0);

comb_tree = cell(size(marker_combinations));
for i_cmb = 1:numel(comb_tree)
	comb_tree{i_cmb} = struct( ...
		'map', {}, ...
		'files', {} ...
		);

	% Allocate structure entry
	comb_tree{i_cmb}(1).map = ...
		containers.Map('KeyType','uint32', 'ValueType','uint32');
	comb_tree{i_cmb}.files = {};
end

% comb_tree contains an array of event times, categorized by combination
% and file.
% comb_tree.map maps the index in data_files to a one-based list.
% comb_tree.files is a cell array indexed by the values of comb_tree.map
% and contains in every field an n-by-2 array of event times in the same order as
% the traces appear in the corresponding files, wherein the column of
% the value for a marker corresponds to the position of the marker in
% marker_combinations.

% Iterate over entries in metadata struct
for i_meta = 1:length(metaS)

	% Get condition token
	curr_cond = metaS(i_meta).condition;

	% Skip if not desired
	if ~isempty(restrict_to_conditions)
		idx_rest_cond = findStrInCstr(restrict_to_conditions, curr_cond, 'i');
		if isempty(idx_rest_cond)
			continue
		end
		if isempty(out_token)
			out_token_temp = sprintf('%s_', restrict_to_conditions{idx_rest_cond});
		else
			out_token_temp = out_token;
		end
	else
		out_token_temp = out_token;
	end

	% Get marker token
	curr_mark = metaS(i_meta).marker;

	% Test whether current entry is relevant
	curr_comb = zeros(1, ncombs);
	for i_comb = 1:ncombs
		curr_comb(i_comb) = ...
			findStrInCstr(marker_combinations{i_comb}, curr_mark, 'i0');
	end

	% Skip if not relevant
	if ~any(curr_comb)
		continue
	end

	% Get relevant parts of file name
	[~, name, ~] = fileparts(dataF{metaS(i_meta).index_f});
	name = textscan(name, '%s', 'Delimiter', '_');
	name = name{1};
	curr_fn = strjoin(name(1:3), '_');
	curr_fn_idx = findStrInCstr(data_files, curr_fn, '0');
	if ~curr_fn_idx
		data_files = [ data_files; curr_fn ];
		curr_fn_idx = length(data_files);
	end

	% Create new comb_traces entry
	for i_comb = find(curr_comb)

		% Register in comb_tree
		keys = comb_tree{i_comb}.map.keys;
		if ~ismember(curr_fn_idx, [keys{:}])
			% Create new file entry in comb_tree{i_comb}.files
			index = length(comb_tree{i_comb}.files) + 1;
			comb_tree{i_comb}.map(curr_fn_idx) = index;
			comb_tree{i_comb}.files{index} = [];
		end

		comb_tree{i_comb}.files{ comb_tree{i_comb}.map(curr_fn_idx) } ...
			(metaS(i_meta).index_t, curr_comb(i_comb)) ...
			= tracesA(metaS(i_meta).index_a, index_tevent) * time_scale;
	end
end

%% Write data in plottable array
corr_data = cell(size(marker_combinations));

for i_cmb = 1:numel(comb_tree)

	% Initialize empty array
	corr_data{i_cmb} = [];

	% Iterate over files
	keys = comb_tree{i_cmb}.map.keys;
	keys = [ keys{:} ];

	for key = keys
		i_file = comb_tree{i_cmb}.map(key);

		% Get array of event times
		time_array = comb_tree{i_cmb}.files{i_file};

		% Skip files where not both columns are populated (~= 0)
		if isempty(time_array) ...
				|| size(time_array, 2) < 2 ...
				|| ~all(any(time_array))
			% No or only one column is populated, skip this file
			continue
		end

		% Detect datasets of unequal lengths
		zero_index = any(time_array' == 0);
		if any(zero_index)
			warning(['Unequal data lengths detected at ' ...
				'combination %s/%s for datafile %s. ' ...
				'Filtering out bad datafile.'], ...
				marker_combinations{i_cmb}{1}, ...
				marker_combinations{i_cmb}{2}, ...
				data_files{key});

			time_array(zero_index,:) = [];	% Delete zero values
		end

		% Append values of this file to corr_data{i_cmb}
		corr_data{i_cmb} = [ corr_data{i_cmb} ; time_array ];

		% DEBUG MESSAGE
		if DEBUG
			fprintf(fid, 'i_comb=% 2d; length=% 5d; file[% 4d]: %s\n', ...
				i_cmb, size(corr_data{i_cmb},1), key, data_files{key});
		end

	end
end

% Close logfile
if DEBUG
	fclose(fid);
end

%% Plot the data
% Define/Initialize values
[len_vert, len_horz] = size(marker_combinations);
n_subplots = len_vert * len_horz;
i_subplot = @(i_row, i_col) i_col + len_horz * (i_row - 1);

color_info = cell(n_subplots, 1);
comb_text = cell(n_subplots, 1);
clusters = cell(n_subplots, 1);
ellipses = cell(n_subplots, 1);
centersOfMass = cell(n_subplots, 1);
timeOrders = zeros([ size(marker_combinations) 3 ]);

% First plot: combinations in separate subplots
fmh = figure;

for i_row = 1:len_vert
	for i_col = 1:len_horz
		i_sub = i_subplot(i_row, i_col);
		ax = subplot(len_vert, len_horz, i_sub);
		this_data = corr_data{i_row, i_col};

		% Get color for scatter plot
		if size(combination_colors) == size(marker_combinations)
			color_info{i_sub} = {'Color', combination_colors{i_row, i_col}};
		elseif numel(combination_colors) >= numel(marker_combinations)
			color_info{i_sub} = {'Color', combination_colors{i_row + len_vert * (i_col - 1)}};
		else
			color_info{i_sub} = [];
		end

		if isempty(this_data)
			% Skip empty combinations (to avoid matrix index error)
			warning('No data found for combination %s/%s.', ...
				marker_combinations{i_row,i_col}{1}, ...
				marker_combinations{i_row,i_col}{2});
			%continue

		else
			% Plot non-empty combination

			% Cluster result, if desired
			clusters{i_sub} = getClusters(this_data, ...
				marker_combinations{i_row, i_col});

			% Calculate error ellipses
			[ellips, com] = getEllipses(this_data, clusters{i_sub});
			for i_el = 1:numel(ellips)
				if any(isnan(ellips{i_el}(:)))
					% Delete ellipse data resulting from too few data points
					ellips{i_el} = [];
					com(i_el,:) = NaN;
				end
			end
			ellipses{i_sub} = ellips;
			centersOfMass{i_sub} = com;

			% Get indices for non-finite event times
			idx_onlyFirst = isfinite(this_data(:,1)) & ~isfinite(this_data(:,2));
			idx_onlySecond = ~isfinite(this_data(:,1)) & isfinite(this_data(:,2));

			% Create a scatter plot
			plot(this_data(:,1), this_data(:,2), '.', color_info{i_sub}{:})
			hold on
			plot(this_data(idx_onlyFirst,1), zeros(sum(idx_onlyFirst),1), '+', color_info{i_sub}{:})
			plot(zeros(sum(idx_onlySecond),1), this_data(idx_onlySecond,2), '+', color_info{i_sub}{:})

			% If clusters have been identified, indicate them
			% (draw noise points as rings)
			p_noise = this_data(~ismember(1:size(this_data,1), vertcat(clusters{i_sub}{:})), :);
			for i_cl = 1:numel(clusters{i_sub})
				%plot(this_data(clusters{i_sub}{i_cl}, 1), ...
				%	this_data(clusters{i_sub}{i_cl}, 2), '.w', 'MarkerSize', 0.5);
				plot(p_noise(:,1), p_noise(:,2), '.w', 'MarkerSize', 0.5);
			end

			% Plot error ellipses
			for i_el = 1:numel(ellips)
				if ~isempty(ellips{i_el})
					% Plot error ellipse
					plot(ellips{i_el}(:,1), ellips{i_el}(:,2), ...
						'-k');%, color_info{i_sub}{:});
					hold on

					% Plot center of mass of error ellipse
					plot(com(i_el,1), com(i_el,2), '+k');
				end
			end
		end

		% Formating and labeling text
		ax.Box = 'On';
% 		ax.YDir = 'reverse';
		ax.XLim = [0 t_max];
		ax.YLim = [0 t_max];

		l = line(ax, xlim, ylim, 'Color', 'k');
		uistack(l, 'bottom')

		xlabel(ax, ['$t_\mathrm{event}(\mathrm{' marker_combinations{i_row,i_col}{1} '})$ [h]'], 'interpreter','latex')
		ylabel(ax, ['$t_\mathrm{event}(\mathrm{' marker_combinations{i_row,i_col}{2} '})$ [h]'], 'interpreter','latex')

		% Get combination description text …
		% (colored output with 'tex' interpreter)
		comb_text{i_sub} = sprintf( ...
			'\\color[rgb]{%f,%f,%f}%s → %s ', ...
			color_info{i_sub}{2}, ...
			marker_combinations{i_row,i_col}{1}, ...
			marker_combinations{i_row,i_col}{2} ...
			);

		% … and add it to the figure
		text(ax, 1, 1, ...
			comb_text{i_sub}, ...
			'Interpreter', 'tex', ...
			'Units', 'normalized', ...
			'HorizontalAlignment', 'right', ...
			'VerticalAlignment', 'top');
	end 
end

% Second plot: all combinations in same axes
fsh = figure;
ax = gca;

for i_row = 1:len_vert
	for i_col = 1:len_horz
		i_sub = i_subplot(i_row, i_col);
		this_data = corr_data{i_row,i_col};

		% Skip empty combinations (to avoid matrix index error)
		if isempty(corr_data{i_row,i_col})
			warning('No data found for combination %s/%s.', ...
				marker_combinations{i_row,i_col}{1}, ...
				marker_combinations{i_row,i_col}{2});
			continue
		end

		% Create a scatter plot
		plot(this_data(:,1), this_data(:,2), '.', color_info{i_sub}{:})
		hold on

		% If clusters have been identified, indicate them
		% (draw noise points as rings)
		p_noise = this_data(~ismember(1:size(this_data,1), vertcat(clusters{i_sub}{:})), :);
		for i_cl = 1:numel(clusters{i_sub})
			plot(p_noise(:,1), p_noise(:,2), '.w', 'MarkerSize', 0.5);
		end

		% Plot error ellipses
		for i_el = 1:numel(ellipses{i_sub})
			if ~isempty(ellipses{i_sub}{i_el})
				plot(ellipses{i_sub}{i_el}(:,1), ellipses{i_sub}{i_el}(:,2), ...
					'-', color_info{i_sub}{:});
			end
		end
	end
end

% Format the figure
% ax.YDir = 'reverse';
ax.XLim = [0 t_max];
ax.YLim = [0 t_max];

l = line(xlim, ylim, 'Color', 'k');
uistack(l, 'bottom')

xlabel('$t_\mathrm{event}^{(n)}$ [h]', 'interpreter','latex')
ylabel('$t_\mathrm{event}^{(n+1)}$ [h]', 'interpreter','latex')

% Add marker description to the figure
text(1, 1, ...
	comb_text, ...
	'Interpreter', 'tex', ...
	'Units', 'normalized', ...
	'HorizontalAlignment','right', ...
	'VerticalAlignment', 'top');


% Third plot: save differences between figures
n_hist_bins = 100;
n_box_bins = 30;
t_x = linspace(0, t_max, 100);

clear('xy_corr')
xy_corr(1:numel(marker_combinations)) = struct();
fphs = cell(1, numel(marker_combinations));
max_t_event = NaN([size(marker_combinations) 2]);

for i_row = 1:len_vert
	for i_col = 1:len_horz
		i_sub = i_subplot(i_row, i_col);
% 		ax = subplot(len_vert, len_horz, i_sub);
		fh_corr = figure;
		ax = axes(fh_corr);
		this_data = corr_data{i_row, i_col};

		% Remove all datapoints that do not have two event times
		idx = any(~isfinite(this_data), 2);
		this_data(idx,:) = [];

		% Prevent indexing error for empty data
		if isempty(this_data)
			continue
		end

		% Get differential data
		data_x = this_data(:,1);
		data_y = this_data(:,2) - this_data(:,1);

		% Get time order for statistics output
		isearlier = data_y > thresh_eq;
		islater = data_y < -thresh_eq;
		timeOrders(i_row, i_col, 1) = ...
			timeOrders(i_row, i_col, 1) + sum(isearlier);
		timeOrders(i_row, i_col, 2) = ...
			timeOrders(i_row, i_col, 2) + sum(~(isearlier | islater));
		timeOrders(i_row, i_col, 3) = ...
			timeOrders(i_row, i_col, 3) + sum(islater);

		% Get histcounts of differential data
		d_data = max(data_x) - min(data_x);
		bin_width_x = d_data / n_hist_bins;
		edges_x = min(data_x) : bin_width_x : max(data_x);
		[counts_x, ~, idx_bin_x] = histcounts(data_x, edges_x);
		area_x = bin_width_x * sum(counts_x);
		dist_x = fitdist(data_x, 'Kernel', 'Kernel', 'normal', 'Support', 'positive');

		d_data = max(data_y) - min(data_y);
		bin_width_y = d_data / n_hist_bins;
		edges_y = min(data_y) : bin_width_y : max(data_y);
		counts_y = histcounts(data_y, edges_y);
		area_y = bin_width_y * sum(bin_width_y);
		dist_y = fitdist(data_y, 'Kernel', 'Kernel', 'normal', 'Support', 'unbound');

		% Get the data indices that belong to each box bin
		box_bin_cat = NaN(length(data_x), 1);
		box_edges = linspace(0, t_max, n_box_bins + 1);
		box_bin_centers = box_edges(1:end-1) + (box_edges(2)-box_edges(1))/2;
		for ic = n_box_bins:-1:1
			idx = find(data_x >= box_edges(ic) & data_x < box_edges(ic+1));
			if ~isempty(idx)
				box_bin_cat(idx) = ic;
			else
				box_bin_centers(ic) = [];
			end
		end

		% Plot differential data
		plot(ax, data_x, data_y, '.', color_info{i_sub}{:});
		hold on

		% Plot distributions
		pdf_mkr1 = pdf(dist_x, t_x) * area_x;
		scale_x = ax.YLim(2) / max(pdf_mkr1) * 0.8;
		plot(ax, t_x, scale_x * pdf_mkr1, '-', color_info{i_sub}{:})

		% Plot boxplot
		boxplot(ax, data_y, box_bin_cat, 'Positions', box_bin_centers, ...
			'Symbol', '', 'PlotStyle', 'traditional', 'Whisker', 0, ...
			'BoxStyle', 'filled', 'MedianStyle', 'target');

		% Format axes
		ax.XLim = [0 t_max];
		ax.Position(3) = ax.Position(3) - 0.1;
		ax_y = axes(fh_corr, 'Position', [ ...
			ax.Position(1) + ax.Position(3), ...
			ax.Position(2), ...
			0.1, ...
			ax.Position(4) ]);
		abscissa_y = linspace(ax.YLim(1), ax.YLim(2), 100);
		pdf_mkr2 = area_y * pdf(dist_y, abscissa_y);
		plot(ax_y, abscissa_y, pdf_mkr2, '-', color_info{i_sub}{:})
		ax_y.YLim = [ 0 1.2*max(pdf_mkr2) ];
% 		ax_y.Box = 'off';
		camroll(ax_y, 90)
		ax_y.XDir = 'reverse';
		ax_y.YAxisLocation = 'origin';
		ax_y.YTick = [];
% 		ax_y.XTickLabel = [];
		ax_y.XAxisLocation = 'top';
% 		ax_y.TickDir = 'both';

		% Plot maxima of probability distributions
		[~, max_ind_x] = max(pdf_mkr1);
		max_t_event(i_row, i_col, 1) = t_x(max_ind_x);
		line(ax, [max_t_event(i_row, i_col, 1) max_t_event(i_row, i_col, 1)], ylim(ax), color_info{i_sub}{:})
		[~, max_ind_y] = max(pdf_mkr2);
		max_t_event(i_row, i_col, 2) = abscissa_y(max_ind_y);
		line(ax_y, [max_t_event(i_row, i_col, 2) max_t_event(i_row, i_col, 2)], ylim(ax_y), color_info{i_sub}{:})
		l = line(ax, ax.XLim, [max_t_event(i_row, i_col, 2) max_t_event(i_row, i_col, 2)], color_info{i_sub}{:});
		uistack(l, 'bottom');

		% Format figure
		ax.XTickMode = 'auto';
		ax.XTickLabelMode = 'auto';
		l = line(ax, ax.XLim, [0 0], 'Color', 'k', 'LineStyle', '-');
		uistack(l, 'bottom')
		xlim(ax_y, ylim(ax))
% 		ax.Box = 'off';
		xlabel(ax, ['$t_\mathrm{event}(\mathrm{' marker_combinations{i_row,i_col}{1} '})$ [h]'], 'interpreter','latex')
		ylabel(ax, ['$t_\mathrm{event}(\mathrm{' marker_combinations{i_row,i_col}{2} '})$ [h]'], 'interpreter','latex')
		title(ax, sprintf('%s → %s', marker_combinations{i_row,i_col}{:}));

		% Copy data for later analysis/debugging
		fphs{i_sub} = fh_corr;

		xy_corr(i_sub).data_x = data_x;
		xy_corr(i_sub).bin_width_x = bin_width_x;
		xy_corr(i_sub).edges_x = edges_x;
		xy_corr(i_sub).counts_x = counts_x;
% 		xy_corr(i_sub).medians_x = medians_x;
		xy_corr(i_sub).box_bin_cat = box_bin_cat;
		xy_corr(i_sub).box_bin_centers = box_bin_centers;
		xy_corr(i_sub).box_edges = box_edges;
		xy_corr(i_sub).area_x = area_x;
		xy_corr(i_sub).dist_x = dist_x;
		xy_corr(i_sub).data_y = data_y;
		xy_corr(i_sub).bin_width_y = bin_width_y;
		xy_corr(i_sub).edges_y = edges_y;
		xy_corr(i_sub).counts_y = counts_y;
		xy_corr(i_sub).area_y = area_y;
		xy_corr(i_sub).dist_y = dist_y;
		xy_corr(i_sub).max_t_event = squeeze(max_t_event(i_row, i_col, :));
	end
end

%% Save figures
% Get time as unique file identifier
date_now = getTime;

% Save first figure (containing several subplots)
save_path = fullfile(out_dir, [ date_now '_' out_token_temp 'eventTimeCorr_multi.pdf' ]);
set(fmh, ...
       'PaperPositionMode', 'manual', ...
       'PaperUnits', 'centimeters', ...
       'PaperSize', [25 15], ...
       'PaperPosition', [0 0 25 15]...
       )
print(fmh, save_path, '-dpdf')

% Save second figure (all combinations within one plot)
save_path = fullfile(out_dir, [ date_now '_' out_token_temp 'eventTimeCorr_single.pdf' ]);
set(fsh, ...
       'PaperPositionMode', 'manual', ...
       'PaperUnits', 'centimeters', ...
       'PaperSize', [15 15], ...
       'PaperPosition', [0 0 15 15]...
       )
print(fsh, save_path, '-dpdf')

% Save third figure (probability distributions)
for i_row = 1:len_vert
	for i_col = 1:len_horz
		save_path = fullfile(out_dir, ...
			sprintf('%s_%seventTimeCorr_%s-%s.pdf', ...
				date_now, out_token_temp, marker_combinations{i_row,i_col}{:}));

		i_sub = i_subplot(i_row, i_col);
		set(fphs{i_sub}, ...
			'PaperPositionMode', 'manual', ...
			'PaperUnits', 'normalized', ...
			'PaperType', 'a4', ...
			'PaperOrientation', 'landscape', ...
			'PaperPosition', [0 0 1 1] ...
			)
		print(fphs{i_sub}, save_path, '-dpdf')
	end
end

%% Finally, print information on how many cells were plotted
% Create statistics matrix
%	statistics_matrix(:,:,1)	number of cells with events for both markers
%	statistics_matrix(:,:,2)	number of cells with event only for first marker
%	statistics_matrix(:,:,3)	number of cells with event only for second marker
%	statistics_matrix(:,:,4)	number of cells without event for any marker
%	statistics_matrix(:,:,5)	total number of cells (sum of the above)
statistics_matrix = zeros([size(marker_combinations) 5], 'uint32');
for i_row = 1:len_vert
	for i_col = 1:len_horz
		is_finite = isfinite(corr_data{i_row,i_col});
		if ~isempty(is_finite)
			statistics_matrix(i_row, i_col, 1) = sum(all(is_finite, 2));
			statistics_matrix(i_row, i_col, 2) = sum( is_finite(:,1) & ~is_finite(:,2) );
			statistics_matrix(i_row, i_col, 3) = sum( ~is_finite(:,1) & is_finite(:,2) );
			statistics_matrix(i_row, i_col, 4) = sum(~any(is_finite, 2));
		else
			statistics_matrix(i_row, i_col, 1:4) = 0;
		end
		statistics_matrix(i_row, i_col, 5) = size(corr_data{i_row,i_col}, 1);
	end
end

% Print statistics
%	“Both”:		found finite event times for both markers for this cell
%	“First”:	found finite event time only for first marker for this cell
%	“Second”:	found finite event time only for second marker for this cell
%	“None”:		found infinite event time/no event for this cell
%	“Total”:	number of all cells in this combination (sum of the above)
logfh = fopen(fullfile(out_dir, [ date_now '_' out_token_temp 'eventTimeCorr.info' ]), 'w');
fprintf(logfh, 'Number of cells plotted:\n');
fprintf(logfh, '%15s  %7s %7s %7s %7s %7s\n', ...
	'Combination', 'Both', 'First', 'Second', 'None', 'Total');
fprintf('Number of cells plotted:\n')
fprintf('%15s  %7s %7s %7s %7s %7s\n', ...
	'Combination', 'Both', 'First', 'Second', 'None', 'Total')
for i_row = 1:len_vert
	for i_col = 1:len_horz
		cmb = marker_combinations{i_row, i_col};
		n_finite = statistics_matrix(i_row, i_col, 1);
		n_finite1 = statistics_matrix(i_row, i_col, 2);
		n_finite2 = statistics_matrix(i_row, i_col, 3);
		n_infinite = statistics_matrix(i_row, i_col, 4);
		n_total = statistics_matrix(i_row, i_col, 5);

		fprintf(logfh, '%7s/%-7s: %7d %7d %7d %7d %7d\n', ...
			cmb{1}, cmb{2}, n_finite, n_finite1, n_finite2, n_infinite, n_total);
		fprintf('%7s/%-7s: %7d %7d %7d %7d %7d\n', ...
			cmb{1}, cmb{2}, n_finite, n_finite1, n_finite2, n_infinite, n_total)
	end
end

% Print information which marker showed signal earlier
%	“First”:	the event indicated by the first marker happens earlier
%	“Equal”:	the events indicated by both markers happen simultaneously
%				(with a tolerance of `2 * thresh_eq`)
%	“Second”:	the event indicated by the second marker happens earlier
fprintf(logfh, '\nEarlier event time:\n');
fprintf(logfh, '%15s  %7s %7s %7s\n', ...
	'Combination', 'First', 'Equal', 'Second');
fprintf('\nEarlier event time:\n');
fprintf('%15s  %7s %7s %7s\n', ...
	'Combination', 'First', 'Equal', 'Second');
for i_row = 1:len_vert
	for i_col = 1:len_horz
		cmb = marker_combinations{i_row, i_col};

		% Normalize `timeOrders`
		timeOrders(i_row, i_col, :) = ...
			timeOrders(i_row, i_col, :) ./ sum(timeOrders(i_row, i_col, :));

		% Print the values
		fprintf(logfh, '%7s/%-7s: %6.1f%% %6.1f%% %6.1f%%\n', ...
			cmb{1}, cmb{2}, 100 * timeOrders(i_row, i_col, :));
		fprintf('%7s/%-7s: %6.1f%% %6.1f%% %6.1f%%\n', ...
			cmb{1}, cmb{2}, 100 * timeOrders(i_row, i_col, :));
	end
end

% Print maxima of probability distributions
fprintf(logfh, '\nMaxima of probability distributions:\n');
fprintf(logfh, '%15s  %7s %7s\n', 'Combination', 'marker1', 'marker2');
fprintf('\nMaxima of probability distributions:\n')
fprintf('%15s  %7s %7s\n', 'Combination', 'marker1', 'marker2')
for i_row = 1:len_vert
	for i_col = 1:len_horz
		cmb = marker_combinations{i_row, i_col};
		fprintf(logfh, '%7s/%-7s: %7.2f %7.2f\n', ...
			cmb{1}, cmb{2}, max_t_event(i_row, i_col, 1), max_t_event(i_row, i_col, 2));
		fprintf('%7s/%-7s: %7.2f %7.2f\n', ...
			cmb{1}, cmb{2}, max_t_event(i_row, i_col, 1), max_t_event(i_row, i_col, 2))
	end
end

% Print centers of mass of error ellipses
fprintf(logfh, '\nCenters of mass of ellipses (x,y):\n');
fprintf('\nCenters of mass of ellipses (x,y):\n');
for i_row = 1:len_vert
	for i_col = 1:len_horz
		i_sub = i_subplot(i_row, i_col);
		com = centersOfMass{i_sub};

		if any(isnan(com))
			% Ignore empty error ellipse
			continue
		end

		cmb = marker_combinations{i_row, i_col};
		for i_el = 1:size(com, 1)
			fprintf(logfh, '%7s/%-7s: (% 3.1f,% 3.1f)\n', ...
				cmb{1}, cmb{2}, com(1), com(2));
			fprintf('%7s/%-7s: (% 3.1f,% 3.1f)\n', ...
				cmb{1}, cmb{2}, com(1), com(2));
		end
	end
end
fclose(logfh);


%% Auxiliary functions
function clusters = getClusters(corr_data, mkr_cmb)
	%% Gets clusters in current dataset
	clusters = {};
	return % Uncomment this line to perform clustering

	if findStrInCstr(mkr_cmb, 'tmrm', 'i0') ...
			&& findStrInCstr(mkr_cmb, 'ros', 'i0')
		cid = dbscan(corr_data, 2.7, 50, 'squaredeuclidean');
		clusters = getLargestClusters(cid, 2);

	elseif findStrInCstr(mkr_cmb, 'lyso', 'i0') ...
			&& findStrInCstr(mkr_cmb, 'tmrm', 'i0')
		cid = dbscan(corr_data, 2.5, 25, 'squaredeuclidean');
		clusters = getLargestClusters(cid, 2);

	elseif findStrInCstr(mkr_cmb, 'lyso', 'i0') ...
			&& findStrInCstr(mkr_cmb, 'ros', 'i0')
		cid = dbscan(corr_data, 2.5, 15, 'squaredeuclidean');
		clusters = getLargestClusters(cid, 1);

	elseif findStrInCstr(mkr_cmb, 'casp', 'i0') ...
			&& findStrInCstr(mkr_cmb, 'toto', 'i0')
		cid = dbscan(corr_data, 2.5, 30, 'squaredeuclidean');
		clusters = getLargestClusters(cid, 1);
	end
end

function clusterPts = getLargestClusters(cid, n)
	%% Gets the n largest clusters in cid
	% Get cluster sizes
	cs = [];
	clu = unique(cid);
	clu = clu(isfinite(clu));
	clu = clu(clu > 0);
	for cli = 1:length(clu)
		cl = clu(cli);
		cs(cl) = length(find(cid == cl));
	end

	% Get n largest clusters
	[~,idx] = sort(cs);
	if length(idx) > n
		idx = idx(1:n);
	end

	% Write indices of points belonging to largest clusters
	clusterPts = cell(n, 1);
	for i_idx = 1:length(idx)
		clusterPts{i_idx} = find(cid == idx(i_idx));
	end
end

function [ellips, centerOfMass] = getEllipses(data, cluster)
	%% Calculates ellipses for data
	% One ellipse per cluster, or if there is no cluster,
	% one ellipse for all datapoints

	nClusters = numel(cluster);

	% Initialize return cell
	if nClusters == 0
		ellips = cell(1);
		centerOfMass = NaN(1, 2);
	else
		ellips = cell(nClusters, 1);
		centerOfMass = NaN(nClusters, 2);
	end

	% Get ellipses
	if nClusters == 0
		[ellips{1}, centerOfMass(1), centerOfMass(2)] = ...
			error_ellipse(data);
		ellips{1} = ellips{1} + centerOfMass;
	else
		for i_el = 1:nClusters
			[ellips{i_el}, centerOfMass(i_el,1), centerOfMass(i_el,2)] = ...
				error_ellipse(data(cluster{i_el},:));
			ellips{i_el} = ellips{i_el} + centerOfMass(i_el,:);
		end
	end
end
