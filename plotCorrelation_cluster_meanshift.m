%plotCorrelation_cluster_meanshift plots correlation scatter plots with
% optional cluster recognition and error ellipses
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
limX = [0 t_max];
limY = [0 t_max];

text_arrow = char(8594);

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

%% Check whether to create error ellipses and whether to cluster
% Check for existing ellipse directive
if exist('makeEllipses', 'var')
	wasEllipseGiven = true;
	if ~islogical(makeEllipses)
		makeEllipses = logical(makeEllipses);
	end
else
	wasEllipseGiven = false;
	makeEllipses = [];
end

% Check for existing clustering directive
if exist('makeClusters', 'var')
	wasClusterGiven = true;
	if ~islogical(makeClusters)
		makeClusters = logical(makeClusters);
	end
else
	wasClusterGiven = false;
	makeClusters = [];
end

% Ask user whether to cluster or show ellipses
if isempty(makeClusters) || isempty(makeEllipses)
	option1 = 'None';
	option2 = 'Only ellipses';
	option3 = 'Ellipses and cluster';
	button = questdlg(['Do you want to show error ellipses of the data?' ...
		'Do you want to cluster the data?'], ...
		'Cluster confirmation', ...
		option1, option2, option3, option3);
	switch button
		case option1
			makeEllipses = false;
			makeClusters = false;
		case option2
			makeEllipses = true;
			makeClusters = false;
		case option3
			makeEllipses = true;
			makeClusters = true;
		otherwise
			return
	end
	clear option1 option2 option3 button
end

%% Prepare ellipse and cluster data
nRows = height(corrT);
clusterT = table(cell(nRows, 1), cell(nRows, 1), cell(nRows, 1), ...
	cell(nRows, 1), cell(nRows, 1), cell(nRows, 1), ...
	'VariableNames', {'center', 'errors', 'labels', 'inEllipse', ...
	'corr_coef', 'corr_stat'});

%% Plot combinations in separate subplots
fhc = cell(1, height(corrT));

for i_cmb = 1:height(corrT)
	% Prepare figure and axes
	fhc{i_cmb} = figure;
	ax = axes(fhc{i_cmb});
	hold(ax, 'on')

	% Formating and labeling axis
	ax.Box = 'On';
	ax.XLim = limX;
	ax.YLim = limY;

	l = line(ax, limX, limY, 'Color', 'k');
	uistack(l, 'bottom')

	xlabel(ax, ['$t_\mathrm{event}(\mathrm{' corrT.markers{i_cmb,1} '})$ [h]'], 'interpreter','latex')
	ylabel(ax, ['$t_\mathrm{event}(\mathrm{' corrT.markers{i_cmb,2} '})$ [h]'], 'interpreter','latex')

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

		% Perform clustering
		clusterT.labels{i_cmb}(corrT.found{i_cmb}(:,1),1) = 1;
		clusterT.labels{i_cmb}(~corrT.found{i_cmb}(:,1),1) = NaN;
		if makeClusters
			% Prepare metadata
			dataName = sprintf('%s: %s %s %s', ...
				strjoin(conditions, '/'), ...
				corrT.markers{i_cmb, 1}, ...
				text_arrow, ...
				corrT.markers{i_cmb, 2});
			fig_names = struct('marker1', corrT.markers{i_cmb, 1}, ...
				'marker2', corrT.markers{i_cmb, 2}, ...
				'cond', conditions);
			if exist('out_dir', 'var')
				fig_names.out_dir = out_dir;
			end

			% Call clustering routine
			[clusterT.labels{i_cmb}(corrT.found{i_cmb}(:,1),1), showBorders] = ...
				doCluster(this_data(corrT.found{i_cmb}(:,1),:), dataName, fig_names);

			if showBorders
				% Get and plot cluster borders
				cluster_borders = voronoi_border(this_data, clusterT.labels{i_cmb});
				plot(ax, cluster_borders(:,1), cluster_borders(:,2), 'k--')
			end
		end

		% Plot error ellipses
		if makeEllipses
			% Get cluster labels
			this_labels = unique(clusterT.labels{i_cmb});
			this_labels = this_labels(isfinite(this_labels));
			nLabels = length(this_labels);

			% Ask user for ellipse asymmetry
			asym_dir = questdlg( ...
				'Which direction of the ellipse shall be asymmetric', ...
				'Ellipse asymmetry', ...
				'Both', 'Large', 'Small', 'Both');
			if isempty(asym_dir)
				asym_dir = 'both';
			end

			% Calculate error ellipse for each cluster
			for iLbl = 1:nLabels
				lbl = this_labels(iLbl);
				idx = ismember(clusterT.labels{i_cmb}, lbl);

				[ellips, comX, comY, ~, ~, ~, ~, isInEllipse] = ...
					error_potato(this_data(idx,:), asym_dir);

				ellips = ellips + [comX comY];
				clusterT.center{i_cmb}(iLbl,:) = [comX, comY];
				clusterT.errors{i_cmb}(iLbl,:) = center_errors(this_data(idx,:));
				clusterT.inEllipse{i_cmb}(iLbl,1) = sum(isInEllipse);

				% Plot error ellipses
				plot(ax, ellips(:,1), ellips(:,2), '-k');
				plot(ax, comX, comY, '+k');

				% Calculate Pearson correlation coefficient
				[corr_coef, n_corr] = pearson_corr(this_data(idx,:));
				clusterT.corr_coef{i_cmb}(iLbl,1) = corr_coef;
				clusterT.corr_stat{i_cmb}(iLbl,1) = n_corr;

				% DEBUG
				% Plot all points inside ellipse as black dots
% 				idx(idx) = isInEllipse;
% 				plot(ax, this_data(idx,1), this_data(idx,2), '.k')
			end
		end
	end

	%% Export figure if requested
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

%% Print information
plotCorrelation_info;

%% Cleanup
if ~wasClusterGiven
	clear makeClusters
end
if ~wasEllipseGiven
	clear makeEllipses
end
clear asym_dir t_max clock_now date_now save_path out_token_temp
clear cluster_borders conditions corr_coef n_corr fig_names
clear i_cmb ax l limX limY showBorders this_data text_dash text_arrow fhc
clear wasClusterGiven wasEllipseGiven comX comY ellips idx dataName
clear nRows iLbl lbl nLabels this_labels isInEllipse

%% Perform clustering, if requested
function [labels, showBorders] = doCluster(data, name, fig_names)
	% Get clustering options (number of clusters and distance measure)
	[bndwdth, kernType, showBorders] = getClusterOpt(name);
	if isempty(bndwdth)
		bndwdth = 1;
	end
	if ~exist('fig_names', 'var')
		fig_names = [];
	end

	% Perform mean shift clustering
	if length(bndwdth) == 1
		% Ignore diagonal for clustering
		labels = runClustering(data, bndwdth);
	else
		% Find separate clusters above and below diagonal
		idxAboveDiag = data(:,1) < data(:,2);
		idxBelowDiag = ~idxAboveDiag;

		labels = NaN(length(idxAboveDiag), 1);
		labels(idxAboveDiag) = runClustering(data(idxAboveDiag,:), bndwdth(2));
		lblOffset = max(unique(labels(isfinite(labels))));
		if isempty(lblOffset)
			lblOffset = 0;
		end
		labels(idxBelowDiag) = lblOffset ...
			+ runClustering(data(idxBelowDiag,:), bndwdth(1));
	end

	% Actual call to clustering function
	function labels = runClustering(D, bndwdth)
		ms = MeanShift(D', bndwdth, kernType);
		labels = ms.cluster();
		[f, ax] = ms.cluster_and_plot(D', true, labels);

		if isstruct(fig_names)
			% Label axes
			cond_str = strjoin(fig_names.cond, '/');
			if ~isempty(cond_str)
				cond_str = [cond_str ': '];
			end
			title(ax, sprintf('%s%s%s%s', cond_str, fig_names.marker1, ...
				char(8594), fig_names.marker2))
			xlabel(ax, sprintf('t_{%s} [h]', fig_names.marker1))
			ylabel(ax, sprintf('t_{%s} [h]', fig_names.marker2))

			% Prepare file name for saving
			cond_str = strjoin(fig_names.cond, '-');
			if ~isempty(cond_str)
				cond_str = strcat(cond_str, char(8211));
			end
			f.UserData.out_dir = fig_names.out_dir;
			f.UserData.out_label = sprintf('%s%s-%s', ...
				cond_str,  fig_names.marker1, fig_names.marker2);
		end
	end

end

%% Calculate Pearson correlation coefficient
function [r_pearson, n_pts] = pearson_corr(points)
%pearson_corr calculates the Pearson correlation coefficient
% Input:
%	points	n-by-2 matrix of points whose correlation to calculate
%
% Returns:
%	r_pearson	scalar Pearson correlation coefficient, in [-1,+1]
%	n_pts		number of points used
%
% Note: If too few (<2) points are given, `r_pearson` is NaN.
%       Non-finite values in `points` may cause undefined behaviour.
	n_pts = size(points, 1);
	if n_pts < 2
		r_pearson = NaN;
		return
	end
	r_pearson = corrcoef(points);
	r_pearson = r_pearson(1,2); % Get off-diagonal element of correlation matrix
end

%% Assess cluster border by Voronoi construction
function cluster_border = voronoi_border(points, labels)
%voronoi_border defines the cluster borders by a voronoi construction

	% Initialize array of cluster border points
	cluster_border = NaN(1,2);

	% Get valid points and cluster information
	idxBad = ~all(isfinite(points), 2);
	points(idxBad,:) = [];
	labels(idxBad) = [];
	[points, ia] = unique(points, 'rows');
	labels = labels(ia);
	labels(isnan(labels)) = 0;
	lbls = unique(labels);
	nLabels = length(lbls);

	% Skip if only one cluster exists
	if nLabels <= 1
		return
	end

	% Skip if too few points for voronoi (to prevent error)
	if size(points, 1) < 3
		return
	end

	% Get indices of cluster members
	idxC = cell(nLabels, 1);
	for iLbl = 1:nLabels
		idxC{iLbl} = find(labels == lbls(iLbl))';
	end

	% Calculate Voronoi pattern
	[vrtx, edge] = voronoin(points);
	border_points = NaN;

	% Identify neighboring points from different clusters
	for i1 = 1:nLabels-1
		idx1 = idxC{i1};
		for i2 = i1+1:nLabels
			idx2 = idxC{i2};

			for j1 = idx1
				for j2 = idx2
					common_border = intersect(edge{j1}, edge{j2}, 'stable');
					if ~isempty(common_border)

						if length(common_border) == 1
							% Unbound cells without common border
							% It must be: `common_border == 1`
							continue

						elseif ismember(1, common_border)
							% Unbound cells with common border
							continue

						elseif ~iscolumn(common_border)
							% Enforce column vector
							common_border = common_border';
						end

						idxB = size(border_points, 1) + (1:length(common_border)+1);
						border_points(idxB, 1) = [common_border; NaN];
					end
				end
			end
		end
	end

	border_points = labelcat(border_points);
	idxB = ~isnan(border_points);
	cluster_border = NaN(length(border_points), 2);
	cluster_border(idxB,:) = vrtx(border_points(idxB),:);
end % end of `voronoi_border`


%% Prompt user for clustering options
function [h, kernType, showBorders] = getClusterOpt(name)
%getClusterOpt assesses clustering options by showing a dialog

	% Process input
	if nargin >= 1
		name = [' in "' name '"'];
	else
		name = '';
	end

	% Build GUI
	win = figure('Name', 'Clustering settings', ...
		'MenuBar', 'none', ...
		'NumberTitle', 'off', ...
		'CloseRequestFcn', @clsWin, ...
		...'WindowStyle', 'modal', ...
		'Units', 'centimeters');
	win.Position(3:4) = [12 7];

	uicontrol(win, ...
		'Style', 'text', ...
		'String', ['Which properties have the clusters' name '?'], ...
		'HorizontalAlignment', 'left', ...
		'Units', 'normalized', ...
		'Position', [.05 .85 .9 .1]);

	% Get first input for number of clusters
	pnl_h1 = uipanel(win, ...
		'Title', 'Bandwidth:', ...
		'BorderType', 'none', ...
		'Units', 'normalized', ...
		'Position', [.3 .7 .4 .15]);
	edt_h1 = uicontrol(pnl_h1, ...
		'Style', 'edit', ...
		'String', '3.5', ...
		'Units', 'normalized', ...
		'Position', [0 0 1 1]);

	% Get second input for number of clusters
	pnl_h2 = uipanel(win, ...
		'Title', 'Bandwidth above diagonal:', ...
		'BorderType', 'none', ...
		'Visible', 'off', ...
		'Units', 'normalized', ...
		'Position', [.05 .7 .4 .15]);
	edt_h2 = uicontrol(pnl_h2, ...
		'Style', 'edit', ...
		'String', '3.5', ...
		'Units', 'normalized', ...
		'Position', [0 0 1 1]);

	% Build radio buttons for diagonal separation
	bg_diag = uibuttongroup(win, ...
		'Title', 'Diagonal treatment', ...
		'SelectionChangedFcn', @toggleEdtNum2, ...
		'BorderType', 'etchedin', ...
		'Units', 'normalized', ...
		'Position', [.05 .4 .425 .25]);
	rb_diag_ign = uicontrol(bg_diag, ...
		'Style', 'radiobutton', ...
		'String', 'Ignore diagonal', ...
		'Units', 'normalized', ...
		'Position', [0 2/3 1 1/3]);
	uicontrol(bg_diag, ...
		'Style', 'radiobutton', ...
		'String', 'Diagonal is border', ...
		'Units', 'normalized', ...
		'Position', [0 1/3 1 1/3]);
	rb_diag_ign.Value = 1;

	% Build radio buttons for distance measure
	clstOpts = ["Epanechnikow", "epanechnikov"; ...
		"Gaussian", "gaussian"];
	nOpts = size(clstOpts, 1);
	bg_kern = uibuttongroup(win, ...
		'Title', 'Kernel type', ...
		'BorderType', 'etchedin', ...
		'Units', 'normalized', ...
		'Position', [.525 .4 .425 .25]);
	for iO = 1:nOpts
		uicontrol(bg_kern, ...
			'Style', 'radiobutton', ...
			'String', clstOpts{iO,1}, ...
			'UserData', clstOpts{iO,2}, ...
			'Max', 1, ...
			'Min', 0, ...
			'Units', 'normalized', ...
			'Position', [0 1-(iO/nOpts) 1 1/nOpts]);
	end
	bg_kern.Children(end).Value = 1;

	chkBorder = uicontrol(win, ...
		'Style', 'checkbox', ...
		'String', 'Show cluster borders', ...
		'Value', true, ...
		'Units', 'normalized', ...
		'Position', [.05 .25 .9 .1] ...
		);

	uicontrol(win, ...
		'Style', 'pushbutton', ...
		'String', 'OK', ...
		'Callback', @clsWin, ...
		'Units', 'normalized', ...
		'Position', [.05 .05 .9 .15]);

	% Wait for user to enter requests
	uiwait(win)

	% Window closing callback
	function clsWin(~, ~, ~)
		h = getBandwidths;
		if rb_diag_ign.Value
			h = h(1);
		end
		showBorders = chkBorder.Value;
		if isempty(bg_kern.SelectedObject)
			kernType = 'epanechnikov';
		else
			kernType = bg_kern.SelectedObject.UserData;
		end
		delete(win);
	end

	% Assess current number of clusters
	function h = getBandwidths
		%getBandwidths reads bandwidths from text entries
		% automatically converts ',' to '.' for better user experience
		h1 = edt_h1.String;
		h1(h1 == ',') = '.';
		h1 = str2double(h1);

		h2 = edt_h2.String;
		h2(h2 == ',') = '.';
		h2 = str2double(h2);

		h = [h1 h2];
	end

	% Show or hide second cluster number input field
	function toggleEdtNum2(~, ~, ~)
		if rb_diag_ign.Value
			% Show only first input field
			pnl_h1.Position(1) = .3;
			pnl_h1.Title = 'Bandwidth:';
			pnl_h2.Visible = 'off';
		else
			% Show both input fields
			pnl_h1.Position(1) = .55;
			pnl_h1.Title = 'Bandwidth below diagonal:';
			pnl_h2.Visible = 'on';
		end
	end

end % end of `getClusterOpt`
