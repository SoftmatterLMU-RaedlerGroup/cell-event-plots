classdef MeanShift
	%MEANSHIFT provides kernel density estimation and mean shift clustering
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

	properties
		S
		h
		K
	end

	methods
		function this = MeanShift(S, h, t)
			%MeanShift constructs a MeanShift instance.
			% Input:
			%	S	sample for spanning the KDE, mxn array,
			%		with number of dimensions m and number of points n
			%	h	bandwidth of kernel
			%	t	(optional) kernel type;
			%		possible values: 'epanechnikov' (default) or 'gaussian'
			this.S = S;
			this.h = h;

			if ~exist('t', 'var')
				this.K = @this.K_epanechnikov;
			else
				switch t
					case 'epanechnikov'
						this.K = @this.K_epanechnikov;
					case 'gaussian'
						this.K = @this.K_gaussian;
					otherwise
						error('Unknown kernel: %s', t)
				end
			end
		end


		function k = kde(this, T)
			%kde calculates a kernel density estimation profile
			% of data this.S at points T with bandwidth this.h
			k = zeros(1, size(T, 2));
			for iT = 1:size(T, 2)
				k(iT) = sum(this.K(sqrt(sum((this.S - T(:,iT)).^2, 1)) ...
					/ this.h));
			end
		end


		function k = kde_mean(this, T)
			%kde calculates a sample mean of data S at points T with bandwidth h
			k = zeros(size(T));
			for iT = 1:size(T, 2)
				k(:,iT) = sum( this.K(sqrt(sum( ...
					(this.S - T(:,iT)).^2, 1)) / this.h) .* this.S, 2 ) ...
					/ this.kde(T(:,iT));
			end
		end


		function [x, y, z] = buildKDE(this, res)
			%buildKDE calculates a KDE grid with given resolution.
			% Input:
			%	res	(optional) grid resolution/number of grid points
			%		If `res` is a scalar, there will be `res` grid points
			%		in both x- and y-direction.
			%		If `res` is a vector of length 2, `res(1)` is the
			%		resolution in x-direction, and `res(2)` in y-direction.
			%		If `res` is empty or not given, a default resolution of
			%		100 for both directions is used.
			%
			% Output:
			%	x	vector of x-values
			%	y	vector of y-values
			%	z	array of KDE values; size: [length(y), length(x)]
			% The output values can be used as an input for `surf`.

			if exist('res', 'var') && ~isempty(res)
				res_x = res(1);
				if length(res) == 2
					res_y = res(2);
				else
					res_y = res(1);
				end
			else
				res_x = 100;
				res_y = 100;
			end

			% Prepare output
			x = linspace(0, 30, res_x);
			y = linspace(0, 30, res_y);
			z = zeros(res_y, res_x);

			% Make density profile
			for ix = 1:res_x
				for iy = 1:res_y
					z(iy,ix) = this.kde([x(ix);y(iy)]);
				end
			end
		end


		function [ax, s] = plotKDE(this, res, ax)
		%plotKDE plots a KDE surface
		% Input:
		%	res	(optional) resolution as required by buildKDE
		%	ax	(optional) axes to plot to; creates new axes if not given
			if ~exist('res', 'var')
				res = [];
			end
	
			% Calculate KDE
			[x, y, z] = this.buildKDE(res);

			% Plot
			if ~exist('ax', 'var')
				ax = axes(figure);
				hold(ax, 'on')
				xlabel(ax, 't1 [h]')
				ylabel(ax, 't2 [h]')
			end
			s = surf(ax, x, y, z, 'FaceColor', 'interp', 'FaceAlpha', 0.75, ...
				'EdgeColor', 'none', 'EdgeAlpha', 0.75);
		end


		function [labels, iterRes] = cluster(this, Told)
			%cluster performs mean-shift clustering
			% Input:
			%	Told		(optional) the points to be clustered, given as
			%				nxm array of m points in n dimensions. If not
			%				given, `S` clustered (blurring).
			% Output:
			%	labels		vector whose i-th element is the label of the
			%				cluster the i-th point from `Told` belongs to.
			%				Labels are contiguously interger numbered,
			%				starting at 1.
			%	iterRes		cell array where the i-th cell contains the
			%				positions of the points after the i-th
			%				iteration

			% Process input
			if ~exist('Told', 'var') || isempty(Told)
				Told = this.S;
			end
			if nargout == 2
				iterRes = cell(1, 1);
				isIterRes = true;
			else
				isIterRes = false;
			end

			% Define thresholds
			iterBreakThresh = 0.02;
			maxDistThresh = 1/6;

			% Perform mean shift
			doLoop = true;
			iRun = 0;
			while doLoop
				iRun = iRun + 1;

				% Calculate mean shift
				T = this.kde_mean(Told);

				% Test for finish
				if all( sqrt(sum((T - Told).^2, 1)) < iterBreakThresh )
					% Iteration finished
					doLoop = false;
				else
					Told = T;
				end

				% Save result of this iteration
				if isIterRes
					iterRes{iRun, 1} = T;
				end
			end

			% Find clusters
			labels = zeros(1, size(T, 2));
			new_lbl = 1;
			while true
				% Find data point without label
				iFree = find(~labels, 1);
				if isempty(iFree)
					break
				end

				% Find similar data points
				iNew = sqrt(sum((T - T(:,iFree)).^2, 1)) <= maxDistThresh;
				lbls = unique(labels(iNew));
				lbls(lbls == 0) = [];

				if isempty(lbls)
					% No label assigned yet, assign new label
					labels(iNew) = new_lbl;
					new_lbl = new_lbl + 1;
				elseif length(lbls) == 1
					% Other label found, use it
					labels(iNew) = lbls;
				else
					% More than one other label found, use smallest
					nl = min(lbls);
					labels(iNew) = nl;

					% Assing new label to all connected clusters
					for l = lbls(lbls ~= nl)
						labels(labels == l) = nl;
					end
				end
			end

			% Reduce cluster number
			lbl_list = unique(labels);
			new_labels = zeros(size(labels));
			for iL = 1:length(lbl_list)
				thisIdx = ismember(labels, lbl_list(iL));
				new_labels(thisIdx) = iL;
			end
			labels = new_labels;
		end


		function [f, ax] = cluster_and_plot(this, T, plot3d, labels)
		%cluster_and_plot performs mean-shift clustering and plots results
		% Input:
		%	T		the points to cluster;
		%			mxn array of n points in m dimensions
		%	plot3d	(optional) Set `true` to plot result 3-dim, else 2-dim
		%	labels	(optional) `labels` as returned from `this.cluster`
		% Output:
		%	f		the figure handle containing the resulting plot
		%	ax		the axis handle containing the resulting plot

			if ~exist('T', 'var') || isempty(T)
				T = this.S;
			end
			if ~exist('plot3d', 'var')
				plot3d = false;
			end
			if ~exist('labels', 'var')
				% Perform mean shift clustering
				labels = this.cluster(T);
			end

			% Build colormap
			clr = label_colors(labels);

			% Plot colored clusters
			f = figure;
			ax = axes(f);
			hold(ax, 'on');

			if plot3d
				% Plot with profile
				[ax, sp] = this.plotKDE([], ax);
				pts = gobjects(1, size(T, 2));
				for iT = 1:size(T, 2)
					pts(iT) = plot3(ax, T(1,iT), T(2,iT), ...
						this.kde(T(:,iT)), '.', 'color', clr(labels(iT), :));
				end

				% Prepare colors for interactive coloring
				f.UserData.colors = struct('col', {}, 'idx', {});
				for iC = 1:size(clr, 1)
					f.UserData.colors(iC).col = clr(iC,:);
					f.UserData.colors(iC).idx = find(ismember(labels, iC));
				end

				% Prepare plot for interactive printing
				% Build context menu for plotting the figure in desired quality
				cm = uicontextmenu;
				submenu = uimenu('Parent', cm, 'Label', 'Surface');
				uimenu('Parent', submenu, 'Label', 'solid', ...
					'Callback', @(~,~,~)set(sp, 'FaceAlpha', 1));
				uimenu('Parent', submenu, 'Label', 'transparent', ...
					'Callback', @(~,~,~)set(sp, 'FaceAlpha', 0.75));
				submenu = uimenu('Parent', cm, 'Label', 'Edges');
				uimenu('Parent', submenu, 'Label', 'show', ...
					'Callback', @(~,~,~)set(sp, 'EdgeColor', 'interp'));
				uimenu('Parent', submenu, 'Label', 'hide', ...
					'Callback', @(~,~,~)set(sp, 'EdgeColor', 'none'));
				uimenu('Parent', submenu, 'Label', 'solid', ...
					'Callback', @(~,~,~)set(sp, 'EdgeAlpha', 1));
				uimenu('Parent', submenu, 'Label', 'transparent', ...
					'Callback', @(~,~,~)set(sp, 'EdgeAlpha', 0.75));
				submenu = uimenu('Parent', cm, 'Label', 'Points');
				uimenu('Parent', submenu, 'Label', 'single-color', ...
					'Callback', {@paint_points, pts, 'mono'});
				uimenu('Parent', submenu, 'Label', 'multi-color', ...
					'Callback', {@paint_points, pts, 'multi'});
				submenu = uimenu('Parent', cm, 'Label', 'Save');
				uimenu('Parent', submenu, 'Label', 'as vector graphics', ...
					'Callback', {@plotThis, 'vector'});
				uimenu('Parent', submenu, 'Label', 'as raster graphics', ...
					'Callback', {@plotThis, 'raster'});
				f.UIContextMenu = cm;

			else
				% Plot flat
				line(ax, [0, 30], [0, 30], 'color', 'k');
				for iT = 1:size(T, 2)
					plot(ax, T(1,iT), T(2, iT), '.', ...
						'color', clr(labels(iT), :))
				end
			end

			% Format axes
			ax.XLim = [0,30];
			ax.YLim = [0,30];
			xlabel(ax, 't1 [h]')
			ylabel(ax, 't2 [h]')
		end

	end


	methods (Static)
		function k = K_gaussian(x)
			%K_gaussian Gaussian kernel
			k = exp(-x.^2);
		end

		function k = K_epanechnikov(x)
			%K_epanechnikov epanechnikov kernel
			k = zeros(size(x));
			idx = abs(x) <= 1;
			k(idx) = 0.75 * (1 - x(idx).^2);
		end
	end

end


function clr = label_colors(labels)
	%label_colors returns a colormap for labelling clusters
	% Input:
	%	labels	vector whose i-th element is the label of the cluster the
	%			i-th point belongs to
	% Output:
	%	clr		colormap with the i-th row being the color for the i-th
	%			cluster. The colormap is based on `lines`. Colors may
	%			repeat. Clusters of size 1 are black.

	% Count labels
	lbl_unique = unique(labels);
	nLbls = length(lbl_unique);
	nLblOccur = zeros(size(lbl_unique));
	for iL = 1:nLbls
		nLblOccur(iL) = sum(ismember(labels, lbl_unique(iL)));
	end

	% Build colormap for clusters
	clr = zeros(nLbls, 3);
	nLblMoreThanOnce = sum(nLblOccur > 1);
	cm = lines(nLblMoreThanOnce);
	iLblMoreThanOnce = 0;
	for iL = 1:nLbls
		if nLblOccur(iL) > 1
			iLblMoreThanOnce = iLblMoreThanOnce + 1;
			clr(iL, :) = cm(iLblMoreThanOnce, :);
		end
	end
end


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
	if ~isfield(gh.UserData, 'out_dir') || ~exist(gh.UserData.out_dir, 'dir')
		gh.UserData.out_dir = '.';
	end

	% Define size of output figure
	gh.Units = 'centimeters';
	gh.PaperPositionMode = 'manual';
	gh.PaperSize = [21 17];
	gh.PaperPosition = [0 0 gh.PaperSize];

	% Assemble file name
	if isfield(gh.UserData, 'out_label')
		save_path = fullfile(gh.UserData.out_dir, [ getTime '_MeanShift_' ...
			gh.UserData.out_label suffix ]);
	else
		[fn, pt] = uiputfile(['*' suffix], 'Save figure name', [gh.UserData.out_dir '/']);
		save_path = fullfile(pt, fn);
	end

	% Save the figure to the given file
	print(gh, renderer{:}, save_path, format)
	fprintf('Figure written to:\n%s\n', save_path);
end


function paint_points(gh, ~, pts, style)
%paint_points paints the points `pts` according to `style`
	if strcmp(style, 'mono')
		for p = pts
			p.Color = 'r';
		end
	elseif strcmp(style, 'multi')
		while ~isa(gh, 'matlab.ui.Figure')
			gh = gh.Parent;
		end
		for iC = 1:numel(gh.UserData.colors)
			for iP = gh.UserData.colors(iC).idx
				pts(iP).Color = gh.UserData.colors(iC).col;
			end
		end
	end
end

