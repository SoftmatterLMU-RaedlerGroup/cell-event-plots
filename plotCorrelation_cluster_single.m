%plotCorrelation_cluster_single plots correlation scatter plots with
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
clusterT = table(cell(nRows, 1), cell(nRows, 1), cell(nRows, 1), cell(nRows, 1), ...
	'VariableNames', {'center', 'errors', 'labels', 'inEllipse'});

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
			dataName = sprintf('%s: %s %s %s', ...
				strjoin(conditions, '/'), ...
				corrT.markers{i_cmb, 1}, ...
				text_arrow, ...
				corrT.markers{i_cmb, 2});
			[clusterT.labels{i_cmb}(corrT.found{i_cmb}(:,1),1), showBorders] = ...
				doCluster(this_data(corrT.found{i_cmb}(:,1),:), dataName);

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

%% Cleanup
if ~wasClusterGiven
	clear makeClusters
end
if ~wasEllipseGiven
	clear makeEllipses
end
clear asym_dir t_max clock_now date_now save_path out_token_temp
clear cluster_borders conditions
clear i_cmb ax l limX limY showBorders this_data text_dash text_arrow fhc
clear wasClusterGiven wasEllipseGiven comX comY ellips idx dataName
clear nRows iLbl lbl nLabels this_labels isInEllipse

%% Perform clustering, if requested
function [labels, showBorders] = doCluster(data, name)
	% Get clustering options (number of clusters and distance measure)
	[num, distMeas, showBorders] = getClusterOpt(name);
	if isempty(num)
		num = 1;
	end

	% Perform clustering
	if length(num) == 1
		% Ignore diagonal for clustering
		labels = runClustering(data, num);
	else
		% Find separate clusters above and below diagonal
		idxAboveDiag = data(:,1) < data(:,2);
		idxBelowDiag = ~idxAboveDiag;

		labels = NaN(length(idxAboveDiag), 1);
		labels(idxAboveDiag) = runClustering(data(idxAboveDiag,:), num(2));
		lblOffset = unique(labels);
		lblOffset = max(lblOffset(isfinite(lblOffset)));
		if isempty(lblOffset)
			lblOffset = 0;
		end
		labels(idxBelowDiag) = lblOffset ...
			+ runClustering(data(idxBelowDiag,:), num(1));
	end

	% Actual call to clustering function
	function labels = runClustering(D, n)
		if n == 1
			labels = ones(size(D,1), 1);
		elseif n == 0
			labels = NaN(size(D,1), 1);
		else
			labels = kmeans(D, n, 'Distance', distMeas, 'Replicates', 10, ...
				'OnlinePhase', 'on', 'EmptyAction', 'drop');
		end
	end

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
function [nClust, clustMeas, showBorders] = getClusterOpt(name)
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
		'String', ['How many clusters are there' name '?'], ...
		'HorizontalAlignment', 'left', ...
		'Units', 'normalized', ...
		'Position', [.05 .85 .9 .1]);

	% Get first input for number of clusters
	pnl_num1 = uipanel(win, ...
		'Title', 'Number of clusters:', ...
		'BorderType', 'none', ...
		'Units', 'normalized', ...
		'Position', [.3 .7 .4 .15]);
	edt_num1 = uicontrol(pnl_num1, ...
		'Style', 'edit', ...
		'String', '1', ...
		'Units', 'normalized', ...
		'Position', [.25 0 .5 1]);
	uicontrol(pnl_num1, ...
		'Style', 'pushbutton', ...
		'String', '+', ...
		'Callback', @(~,~)chNCl(edt_num1, 1), ...
		'Units', 'normalized', ...
		'Position', [.75 0 .25 1]);
	uicontrol(pnl_num1, ...
		'Style', 'pushbutton', ...
		'String', '-', ...
		'Callback', @(~,~)chNCl(edt_num1, -1), ...
		'Units', 'normalized', ...
		'Position', [0 0 .25 1]);

	% Get second input for number of clusters
	pnl_num2 = uipanel(win, ...
		'Title', 'Clusters above diagonal:', ...
		'BorderType', 'none', ...
		'Visible', 'off', ...
		'Units', 'normalized', ...
		'Position', [.05 .7 .4 .15]);
	edt_num2 = uicontrol(pnl_num2, ...
		'Style', 'edit', ...
		'String', '0', ...
		'Units', 'normalized', ...
		'Position', [.25 0 .5 1]);
	uicontrol(pnl_num2, ...
		'Style', 'pushbutton', ...
		'String', '+', ...
		'Callback', @(~,~)chNCl(edt_num2, 1), ...
		'Units', 'normalized', ...
		'Position', [.75 0 .25 1]);
	uicontrol(pnl_num2, ...
		'Style', 'pushbutton', ...
		'String', '-', ...
		'Callback', @(~,~)chNCl(edt_num2, -1), ...
		'Units', 'normalized', ...
		'Position', [0 0 .25 1]);

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
	clstOpts = ["Squared euclidean", "sqeuclidean"; ...
		"Cosine", "cosine"; ...
		"Correlation", "correlation"];
	nOpts = size(clstOpts, 1);
	bg_dist = uibuttongroup(win, ...
		'Title', 'Distance measure', ...
		'BorderType', 'etchedin', ...
		'Units', 'normalized', ...
		'Position', [.525 .4 .425 .25]);
	for iO = 1:nOpts
		uicontrol(bg_dist, ...
			'Style', 'radiobutton', ...
			'String', clstOpts{iO,1}, ...
			'UserData', clstOpts{iO,2}, ...
			'Max', 1, ...
			'Min', 0, ...
			'Units', 'normalized', ...
			'Position', [0 1-(iO/nOpts) 1 1/nOpts]);
	end
	bg_dist.Children(end).Value = 1;

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
		nClust = getNClusters;
		if rb_diag_ign.Value
			nClust = nClust(1);
		end
		showBorders = chkBorder.Value;
		if isempty(bg_dist.SelectedObject)
			clustMeas = 'cosine';
		else
			clustMeas = bg_dist.SelectedObject.UserData;
		end
		delete(win);
	end

	% Assess current number of clusters
	function nCl = getNClusters
		nCl1 = uint8(str2double(edt_num1.String));
		nCl2 = uint8(str2double(edt_num2.String));
		nCl = [nCl1 nCl2];
	end

	% Change number of clusters
	function chNCl(edt_num, delta)
		nCl = getNClusters;
		if edt_num == edt_num1
			nCl = nCl(1);
		else
			nCl = nCl(2);
		end
		nCl = nCl + delta;
		edt_num.String = num2str(nCl);
	end

	% Show or hide second cluster number input field
	function toggleEdtNum2(~, ~, ~)
		nCl = double(getNClusters);
		if rb_diag_ign.Value
			% Show only first input field
			pnl_num1.Position(1) = .3;
			pnl_num1.Title = 'Number of clusters:';
			pnl_num2.Visible = 'off';
			chNCl(edt_num1, sum(nCl) - nCl(1))
		else
			% Show both input fields
			pnl_num1.Position(1) = .55;
			pnl_num1.Title = 'Clusters below diagonal:';
			pnl_num2.Visible = 'on';

			delta = nCl(1) / 2;
			chNCl(edt_num1, ceil(delta) - nCl(1))
			chNCl(edt_num2, floor(delta) - nCl(2))
		end
	end

end % end of `getClusterOpt`
