%%
% If the workspace holds a logical variable `use_Aquad`:
%	true:	use a_quad for ROS
%	false:	use initial slope for ROS [default]
%
% If the workspace holds a logical variable `invert_slope`:
%	true:	use reciprocal slope [default]
%	false:	use direct value of slope
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

%% Set constants and parameters
IDX_A_QUAD = 2;
IDX_SLOPE = 7;

if ~exist('use_Aquad', 'var')
	use_Aquad = false;
else
	use_Aquad = logical(use_Aquad);
end

if ~exist('invert_slope', 'var')
	invert_slope = false;
else
	invert_slope = logical(invert_slope);
end

if ~exist('makeEllipse', 'var')
	makeEllipse = false;
else
	makeEllipse = logical(makeEllipse);
end

if ~exist('manualAxLim', 'var')
	manualAxLim = true;
else
	manualAxLim = logical(manualAxLim);
end

[use_Aquad, invert_slope, xlog, ylog, makeEllipse, fit_method, manualAxLim] = ...
	getAxisScales(use_Aquad, invert_slope, makeEllipse, manualAxLim);

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
	save_path = fullfile(out_dir, [getTime '_' out_token_temp 'eventCorr_ROS']);
else
	save_path = false;
end

%% Get ROS row indices
rowsRos = find(any(ismember(corrT.markers, 'ros'), 2));
nRosRows = length(rowsRos);
fhc = cell(1, nRosRows);
info_data = table(strings(nRosRows, 1), NaN(nRosRows, 2), NaN(nRosRows, 1), zeros(nRosRows, 1), ...
		'VariableNames', {'markers', 'center_of_mass', 'corr_coef', 'number'});
text_arrow = char(8594);

%% Plot and export
for iRos = 1:nRosRows
	rosRow = rowsRos(iRos);
	rosCol = ismember(corrT.markers(rosRow,:), 'ros');
	otherCol = find(~rosCol);
	rosCol = find(rosCol);

	%% Prepare figure and axes
	fhc{iRos} = figure;
	ax = axes(fhc{iRos});
	hold(ax, 'on')

	% Formating and labeling axis
	ax.Box = 'On';
%	ax.XLim = [0 t_max];
	ax.YLim = [0 t_max];

	if xlog
		ax.XScale = 'log';
	else
		ax.XScale = 'linear';
	end

	if ylog
		ax.YScale = 'log';
	else
		ax.YScale = 'linear';
	end

	if use_Aquad
		x_quantity = 'a_\mathrm{quad}';
	else
		x_quantity = '\alpha_0';
	end

	if invert_slope
		x_quantity = strcat(x_quantity, '^{-1}');
	end

	xlabel(ax, ['$' x_quantity '(\mathrm{' corrT.markers{rosRow,rosCol} '})$ [a.u.]'], 'interpreter','latex')
	ylabel(ax, ['$t_\mathrm{event}(\mathrm{' corrT.markers{rosRow,otherCol} '})$ [h]'], 'interpreter','latex')

	title(ax, ...
		sprintf( '\\color[rgb]{%f,%f,%f}%s: %s %s %s', ...
			corrT.color(rosRow,:), ...
			strjoin(conditions, '/'), ...
			corrT.markers{rosRow, rosCol}, ...
			text_arrow, ...
			corrT.markers{rosRow, otherCol}) ...
		)

	%% Plot the event data
	% Get times and slopes to plot
	thisGoodTimes = corrT.found{rosRow}(:,1);
	thisBadTimes = ~(thisGoodTimes | corrT.found{rosRow}(:,1+otherCol));
	this_times = corrT.t_event{rosRow}(thisGoodTimes,otherCol);

	if use_Aquad
		this_slopes = corrT.params{rosRow,rosCol}(:,IDX_A_QUAD);
	else
		this_slopes = corrT.state{rosRow,rosCol}(:,IDX_SLOPE);
	end

	% Invert slope
	if invert_slope
		this_slopes = 1 ./ this_slopes;
	end

	% Fit guide for the eye into data
	maxXgood = max(this_slopes(thisGoodTimes));
	switch(fit_method)
		case 'exp'
			% Fit exponential decay into data
			expdec = @(x, A, B, C) A * exp(-x / B) + C;
			p = fit_mle(expdec, this_slopes(thisGoodTimes), this_times, ...
				'InitialParameters', [30, 3, 2], ...
				'ParameterBounds', [0, 0.01, -30; 60, 50, 30], ...
				'ObjectiveFcn', @objfun);
			%fprintf('A=%f, B=%f, C=%f\n', p)
			x_fit = linspace(0, maxXgood, 20 * maxXgood);
			y_fit = expdec(x_fit, p(1), p(2), p(3));

		case 'hyp'
			% Fit hyperbolic descent
			hypdesc = @(x, A) A * x.^(-1);
			p = fit_mle(hypdesc, this_slopes(thisGoodTimes), this_times, ...
				'InitialParameters', 20, ...
				'ParameterBounds', [0; 60], ...
				'ObjectiveFcn', @objfun);
			x_fit = linspace(0, maxXgood, 20 * maxXgood);
			y_fit = hypdesc(x_fit, p(1));

		case 'const'
			% Fit constant level
			constlev = @(~,A) A;
			p = fit_mle(constlev, this_slopes(thisGoodTimes), this_times, ...
				'InitialParameters', mean(this_times), ...
				'ParameterBounds', [min(this_times); max(this_times)], ...
				'ObjectiveFcn', @(~, ~, dep, A, ~) sum(abs(dep - A).^2));
			x_fit = linspace(0, maxXgood, 20 * maxXgood);
			y_fit = x_fit;
			y_fit(:) = p;

		otherwise
			% No fitting
			x_fit = [];
			y_fit = [];
	end

	% Create and plot error ellipse
	if makeEllipse && size(this_times, 1) > 1
		elps_times = this_times;
		elps_slopes = this_slopes(thisGoodTimes);

		if xlog
			elps_slopes = log10(elps_slopes);
		end
		if ylog
			elps_times = log10(elps_times);
		end

		[ellips, comX, comY] = error_ellipse([elps_slopes elps_times]);
		ellips = ellips + [comX comY];

		if xlog
			ellips(:,1) = 10 .^ ellips(:,1);
			comX = 10 .^ comX;
		end
		if ylog
			ellips(:,2) = 10 .^ ellips(:,2);
			comY = 10 .^ comY;
		end

	else
		ellips = zeros(0, 2);
		comX = NaN;
		comY = NaN;
	end

	% Plot the data
	plot(ax, ...
		this_slopes(thisGoodTimes), ...
		this_times, ...
		'.', 'Color', corrT.color(rosRow,:))
	plot(ax, ...
		x_fit, y_fit, '-r')
	plot(ax, ...
		this_slopes(thisBadTimes), ...
		repmat(ax.YLim(1), sum(thisBadTimes), 1), ...
		'+', 'Color', corrT.color(rosRow,:))
	plot(ax, ellips(:,1), ellips(:,2), '-k')
	plot(ax, comX, comY, '+k')

	% Set axis limits for logarithmic plot
	maxS = max(this_slopes(thisGoodTimes));
	if xlog
		if ax.XLim(2) < maxS
			ax.XLim(2) = maxS;
		end

		minS = min(this_slopes(this_slopes > 0 & thisGoodTimes));
		if ax.XLim(1) > minS
			ax.XLim(1) = minS;
		end
	else
		if ax.XLim(2) > maxS
			ax.XLim(2) = maxS;
		end
		if ax.XLim(1) < 0
			ax.XLim(1) = 0;
		end
	end
	if ylog
		minT = min(this_times);
		if ax.YLim(1) < minT
			ax.YLim(1) = minT;
		end
	else
		if ax.YLim(1) ~= 0
			ax.YLim(1) = 0;
		end
		if ax.YLim(2) ~= 30
			ax.YLim(2) = 30;
		end
	end

	%% Prompt user for axis limits, if requested
	if manualAxLim
		newLim = inputdlg( ...
			{'Lower x axis limit:', 'Upper x axis limit:', ...
				'Lower y axis limit:', 'Upper y axis limit:'}, ...
			sprintf('Choose axis limits for %s%s %s %s%s', ...
				char(8220), corrT.markers{rosRow, 1}, text_arrow, ...
				corrT.markers{rosRow, 2}, char(8221)), ...
			1, ...
			cellstr(string([ax.XLim, ax.YLim])), ...
			struct('Resize', 'on', 'WindowStyle', 'normal'));

		% Apply user requests
		if ~isempty(newLim)
			% Lower x axis limit
			newXll = str2double(strrep(newLim{1}, ',', '.'));
			if newXll == ax.XLim(1)
				% do nothing
			elseif isfinite(newXll) && newXll > 0
				ax.XLim(1) = newXll;
			else
				warning('Ignoring invalid lower x axis limit: %s', newLim{1});
			end

			% Upper x axis limit
			newXlu = str2double(strrep(newLim{2}, ',', '.'));
			if newXlu == ax.XLim(2)
				% do nothing
			elseif isfinite(newXlu) && newXlu > 0
				ax.XLim(2) = newXlu;
			else
				warning('Ignoring invalid upper x axis limit: %s', newLim{2});
			end

			% Lower y axis limit
			newYll = str2double(strrep(newLim{3}, ',', '.'));
			if newYll == ax.YLim(1)
				% do nothing
			elseif isfinite(newYll) && newYll > 0
				ax.YLim(1) = newYll;
			else
				warning('Ignoring invalid lower y axis limit: %s', newLim{3});
			end

			% Upper y axis limit
			newYlu = str2double(strrep(newLim{4}, ',', '.'));
			if newYlu == ax.YLim(2)
				% do nothing
			elseif isfinite(newYlu) && newYlu > 0
				ax.YLim(2) = newYlu;
			else
				warning('Ignoring invalid upper y axis limit: %s', newLim{4});
			end
		end
	end

	% Calculate Pearson correlation coefficient
	r_pearson = corrcoef(this_slopes(thisGoodTimes), this_times);
	r_pearson = r_pearson(1,2); % Get off-diagonal element of correlation matrix

	%% Export figure to file
	if ischar(save_path)
		% Export plot
		save_path_plot = [save_path '.ps'];
		set(fhc{iRos}, ...
			'PaperPositionMode', 'manual', ...
			'PaperUnits', 'centimeters', ...
			'PaperSize', [17 15], ...
			'PaperPosition', [0 0 17 15] ...
			)
		print(fhc{iRos}, '-painters', save_path_plot, '-dpsc', '-append')

		% Export info
		info_data.markers(iRos) = sprintf('%s %s %s', ...
			corrT.markers{rosRow, rosCol}, text_arrow, corrT.markers{rosRow, otherCol});
		info_data.center_of_mass(iRos,:) = [comX, comY];
		info_data.corr_coef(iRos) = r_pearson;
		info_data.number(iRos) = sum(thisGoodTimes);
	end
end

%% Export data information
if ischar(save_path)
	save_path_info = [save_path '.info'];
	logfh = fopen(save_path_info, 'w');

	% Print correlation coefficients
	fprintf(logfh, 'Pearson correlation coefficient');
	fprintf('Pearson correlation coefficient');
	fprintf(logfh, '\n%15s %6s %11s\n', 'Combination', 'Number', 'Correlation');
	fprintf('\n%15s %6s %11s\n', 'Combination', 'Number', 'Correlation');

	for iRos = 1:nRosRows
		fprintf(logfh, '%14s: %6d %11.2f\n', ...
			info_data{iRos,'markers'}, info_data{iRos,'number'}, info_data{iRos,'corr_coef'});
		fprintf('%14s: %6d %11.2f\n', ...
			info_data{iRos,'markers'}, info_data{iRos,'number'}, info_data{iRos,'corr_coef'});
	end

	% (Optional) Print ellipse centers
	if makeEllipse
		fprintf(logfh, '\nCenters of error ellipses');
		fprintf('\nCenters of error ellipses');
		fprintf(logfh, '\n%15s %6s %6s\n', 'Combination', 'Slope', 'Time');
		fprintf('\n%15s %6s %6s\n', 'Combination', 'Slope', 'Time');

		for iRos = 1:nRosRows
			fprintf(logfh, '%14s: %6.2f %6.2f\n', ...
				info_data{iRos,'markers'}, info_data{iRos,'center_of_mass'});
			fprintf('%14s: %6.2f %6.2f\n', ...
				info_data{iRos,'markers'}, info_data{iRos,'center_of_mass'});
		end
	end

	fclose(logfh);
end

%% Cleanup
clear ax comX comY conditions ellips elps_slopes elps_times
clear IDX_A_QUAD IDX_SLOPE info_data iRos logfh maxS minS minT nRosRows
clear otherCol out_token_temp rosCol rosRow rowsRos r_pearson save_path
clear save_path_info save_path_plot t_max text_arrow this_slopes this_times
clear thisBadTimes thisGoodTimes xlog ylog y_quantity

% end of main part of `plotCorrelation_ROS`

%% Auxiliary function
function [aquad, inv, xlg, ylg, mkEl, ftm, mAxLm] = ...
	getAxisScales(aquad, inv, mkEl, mAxLm)
%getAxisScales assesses the axis scales by showing a dialog
%
% Meanings of parameters (logical):
%	aquad	Plot quadratic parabola coefficient instead of tangent slope
%	inv		Plot inverse slope
%	xlg		Scale x-axis logarithmically instead of linearly
%	ylg		Scale y-axis logarithmically instead of linearly
%	mkEl	Calculate and plot error ellipse
%	ftm		Fit method; one of: 'none', 'exp', 'hyp', 'const'
%	mAxLm	Prompt user for axis limits

	% Get default values
	xlg = false;
	ylg = false;
	ftm = 'none';

	if aquad
		ylg = true;
		xlg = true;
	end

	if inv
		xlg = true;
		ylg = true;
	end

	% Build GUI
	win = figure('Name', 'Plot properties', ...
		'MenuBar', 'none', ...
		'NumberTitle', 'off', ...
		'CloseRequestFcn', @clsWin, ...
		'WindowStyle', 'modal', ...
		'Units', 'centimeters');
	win.Position(3:4) = [14 5];

	pnl = uipanel(win, ...
		'Title', 'Plot settings:', ...
		'Units', 'normalized', ...
		'Position', [.02 .52 .51 .46]);

	btn_inv = uicontrol(pnl, ...
		'Style', 'checkbox', ...
		'String', 'Invert slope', ...
		'Value', inv, ...
		'Units', 'normalized', ...
		'Position', [.02 .75 .96 .25]);

	btn_aquad = uicontrol(pnl, ...
		'Style', 'checkbox', ...
		'String', 'Use quadratic parabola coefficient', ...
		'Value', aquad, ...
		'Units', 'normalized', ...
		'Position', [.02 .5 .96 .25]);

	btn_elps = uicontrol(pnl, ...
		'Style', 'checkbox', ...
		'String', 'Make error ellipses', ...
		'Value', mkEl, ...
		'Units', 'normalized', ...
		'Position', [.02 .25 .96 .25]);

	btn_mAxLm = uicontrol(pnl, ...
		'Style', 'checkbox', ...
		'String', 'Set axis limits manually', ...
		'Value', mAxLm, ...
		'Units', 'normalized', ...
		'Position', [.02 0 .96 .25]);

	ftbtg = uibuttongroup(win, ...
		'Title', 'Fit settings:', ...
		'Units', 'normalized', ...
		'Position', [.57 .52 .41 .46]);

	uicontrol(ftbtg, ...
		'Style', 'radiobutton', ...
		'String', 'None', ...
		'Value', true, ...
		'UserData', 'none', ...
		'Units', 'normalized', ...
		'Position', [.02 .75 .96 .25]);

	uicontrol(ftbtg, ...
		'Style', 'radiobutton', ...
		'String', 'Exponential decay', ...
		'Value', false, ...
		'UserData', 'exp', ...
		'Units', 'normalized', ...
		'Position', [.02 .5 .96 .25]);

	uicontrol(ftbtg, ...
		'Style', 'radiobutton', ...
		'String', 'Hyperbolic descent', ...
		'Value', false, ...
		'UserData', 'hyp', ...
		'Units', 'normalized', ...
		'Position', [.02 .25 .96 .25]);

	uicontrol(ftbtg, ...
		'Style', 'radiobutton', ...
		'String', 'Constant level', ...
		'Value', false, ...
		'UserData', 'const', ...
		'Units', 'normalized', ...
		'Position', [.02 0 .96 .25]);

	pnl = uipanel(win, ...
		'Title', 'Logarithmic axes:', ...
		'Units', 'normalized', ...
		'Position', [.02 .27 .96 .21]);

	btn_ylg = uicontrol(pnl, ...
		'Style', 'checkbox', ...
		'String', 'Logarithmic y-Axis', ...
		'Value', ylg, ...
		'Units', 'normalized', ...
		'Position', [0.02 .1 .46 .8]);

	btn_xlg = uicontrol(pnl, ...
		'Style', 'checkbox', ...
		'String', 'Logarithmic x-Axis', ...
		'Value', xlg, ...
		'Units', 'normalized', ...
		'Position', [0.52 .1 .46 .8]);

	uicontrol(win, ...
		'Style', 'pushbutton', ...
		'String', 'OK', ...
		'Callback', @clsWin, ...
		'Units', 'normalized', ...
		'Position', [.02 .02 .96 .21]);

	% Wait for user to enter preferences
	uiwait(win)

	% Window close function
	function clsWin(~, ~, ~)
		inv = btn_inv.Value;
		aquad = btn_aquad.Value;
		xlg = btn_xlg.Value;
		ylg = btn_ylg.Value;
		mkEl = btn_elps.Value;
		ftm = ftbtg.SelectedObject.UserData;
		mAxLm = btn_mAxLm.Value;

		delete(win)
	end
end

function dist_sum = objfun(fun, indep, dep, ~, ~)
%objfun is the objective function for fitting the guide-for-the-eye exp
	n_eval = 500;
	n_data = length(indep);

	% Calculate values of model function
	x_eval = linspace(min(indep), max(indep), n_eval);
	y_eval = fun(x_eval);

	% Iterate through points and calculate distance
	dist_sum = 0;
	for i = 1:n_data
		dist_sum = dist_sum + ...
			min(sqrt( (indep(i) - x_eval).^2 + (dep(i) - y_eval).^2 ));
	end
end