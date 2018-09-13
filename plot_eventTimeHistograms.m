% plot_eventTimeHistograms.m
% This script plots histograms of event times for the different markers.
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
markers = {'lyso', 'ros', 'tmrm', 'casp', 'psiva', 'pi'};%, 'toto'};
marker_names = {'LysoTracker', 'CellRox', 'TMRM', ...
	'Cell Event Caspase 3/7', 'pSIVA-IANBD', 'Nucleus (PI)'};%, ...
% 	'Nucleus (toto)'};
bar_colors = [0 0.8 0; 0.8 1 0; 1 1 0; 1 0.6 0; 1 0 0; 0.4 0 0; 0.4 0 0];

index_eventtime = 2;		% Index of event time column
time_scale = 1;				% time scaling factor (for conversion to hours)
time_scale0 = 1 / 6;		% time scaling factor (for conversion from frames to hours)
n_bins = 40;					% Number of bins in histogram plots
y_label_hist = 'Cell count';	% Label of y-axis in histogram plots
x_label_hist = 'Time [h]';		% Label of x-axis in histogram plots

% Define white list of conditions to be used
if ~exist('restrict_to_conditions', 'var')
	restrict_to_conditions = {};
end

% Define whether to append new markers to the list
if ~exist('appendNewMarkers', 'var') || ~islogical(appendNewMarkers)
	appendNewMarkers = false;
end

% Define additional token to be added to output file names
if ~exist('out_token', 'var')
	out_token = [];
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

%% Prepare axes
% Get markers and event times in current measurement
if ~exist('markers', 'var') || ~iscell(markers')
	markers = cell(0);
end
if ~exist('marker_names', 'var') || numel(marker_names) ~= numel(markers)
	marker_names = markers;
end
event_times = cell(0);
time_vec_max = cell(0);

for i_trace = 1:length(metaS)

	% Get condition token
	curr_cond = metaS(i_trace).condition;

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

	% Assess marker of current trace
	current_marker = metaS(i_trace).marker;
	i_marker = findStrInCstr(markers, current_marker, 'i');

	if isempty(i_marker)
		if appendNewMarkers
			% If marker not registered yet, append it
			i_marker = length(markers) + 1;
			markers{i_marker} = current_marker;
		else
			continue
		end
	end

	% Write event time of trace
	if numel(event_times) < i_marker
		event_times{i_marker} = [];
	end
	event_times{i_marker} = ...
		[ event_times{i_marker} tracesA(i_trace, index_eventtime) ];

	% Write number of frames of trace
	if numel(time_vec_max) < i_marker
		time_vec_max{i_marker} = [];
	end
	time_vec_max{i_marker} = ...
		[ time_vec_max{i_marker} get_nFrames(metaS(i_trace).measurement, time_scale0) ];
end

% Skip if no cells were found
if isempty(event_times) || isempty(time_vec_max)
	return
end

% Ensure consistent length of event_times and time_vec_max
if numel(event_times) < numel(markers) || numel(time_vec_max) < numel(markers)
	for i_mkr = 1:numel(markers)
		if i_mkr > numel(event_times)
			event_times{i_mkr} = [];
		end
		if i_mkr > numel(time_vec_max)
			time_vec_max{i_mkr} = [];
		end
	end
end

% Perform time scaling on time-dependent data
for i_marker = 1:numel(markers)
	event_times{i_marker} = time_scale * event_times{i_marker};
	time_vec_max{i_marker} = time_scale * time_vec_max{i_marker};
end

% Initialize figure with axes
fh = figure;
aha = tight_subplot(numel(markers), 1, 0, .1, .1);
hha = cell(0);
bar_center = @(edges, i_bin) (edges(i_bin+1) + edges(i_bin)) / 2;

% Initialize fit options
% gauss = @(A, mu, sig, x) A * exp( - (x-mu).^2 /2 /sig^2 );
% gamma_dist = @(p, b, x) b.^p / gamma(p) * (x.^(p-1) .* exp(-b.*x));
% burr_dist = @(c, k, x) c * k * ( x.^(c-1) ) ./ ( (1+(x.^c)) .^ (k+1) );
logN = @(a, mu, sigma, x) a * ( exp(- (log(x)-mu).^2 /2 /sigma^2 ) ./x ) /sigma /sqrt(2*pi);
logNormPdf = fittype(logN);
logNormOpt = fitoptions(logNormPdf);
logNormOpt = fitoptions(logNormOpt, 'Lower', [0 eps eps], 'StartPoint', [50 15 10]);
fit_result = cell(0);

%% Plot
% Get limit for X axis
event_time_max = double(max([time_vec_max{:}]));
pdf_max_t = NaN(numel(markers), 1);
hist_data = struct([]);

for i_marker = 1:numel(markers)
	% Plot histogram
	hha{i_marker} = histogram(aha(i_marker), event_times{i_marker}, 0 : event_time_max / n_bins : event_time_max);
	hha{i_marker}.EdgeColor = [ 1 1 1 ];
	hha{i_marker}.EdgeAlpha = 1;
	hha{i_marker}.FaceAlpha = 1;
	hha{i_marker}.LineWidth = 1;
	hha{i_marker}.BinLimits = [0, event_time_max];
% 	hha{i_marker}.BinEdges = 0 : event_time_max / n_bins : event_time_max;
	if exist('bar_colors', 'var') && size(bar_colors, 1) >= i_marker && size(bar_colors, 2) == 3
		hha{i_marker}.FaceColor = bar_colors(i_marker,:);
	end
	hold(aha(i_marker), 'on');

	% Adjust axes
	cell_count_max = max(hha{i_marker}.Values);
	x_max = 1.1 * event_time_max;
	y_max = 1.2 * cell_count_max;
	if x_max == 0
		x_max = 1;
	end
	if y_max == 0
		y_max = 1;
	end
	axis(aha(i_marker), [ 0 x_max 0 y_max ]);

	set(aha(i_marker), ...
		'Box', 'on', ...
		'XTickLabelMode', 'auto', ...
		'YTickLabelMode', 'auto', ...
		'XLimMode', 'manual', ...
		'YLimMode', 'manual');
	ylabel(aha(i_marker), y_label_hist);
	text(aha(i_marker), ...
		0, 1, ...
		[ ' ' marker_names{i_marker} ], ...
		'Units', 'normalized', ...
		'HorizontalAlignment','left', ...
		'VerticalAlignment', 'top');

	if i_marker == length(aha)
		% Lowest plot
		xlabel(aha(i_marker), x_label_hist);
	else
		% Above lowest plot
		set(aha(i_marker), 'XTickLabels', '');
	end

	% Fit data if not empty
	fit_y = hha{i_marker}.Values(:);
	fit_x = bar_center(hha{i_marker}.BinEdges, 1:length(fit_y));
	fit_x = fit_x(:);
	area = hha{i_marker}.BinWidth * sum(fit_y);

	if ~isempty(fit_x) && ~isempty(time_vec_max{i_marker})
		% Prepare fit
		fit_upper_time = max(time_vec_max{i_marker});

		% Fit gauss to histogram
% 		fit_result{i_marker} = fit(fit_x, fit_y, ft, fo);
% 		pdf_mkr = fit_result{i_marker}(fit_x) * area;

		% Fit lognormal distribution to histogram
		fit_result{i_marker} = fit(fit_x, fit_y, logNormPdf, logNormOpt);
		pdf_mkr = fit_result{i_marker}(fit_x);... * area;

		% Use kernel density estimator
% 		fit_data = hha{i_marker}.Data(:);
% 		fit_data = fit_data(isfinite(fit_data));
% 		fit_result{i_marker} = fitdist(fit_data, 'Kernel', 'Kernel', 'epanechnikov', 'Support', 'positive');
% 		pdf_mkr = pdf(fit_result{i_marker}, fit_x) * area;
		
		[~, max_idx_pdf] = max(pdf_mkr);
		pdf_max_t(i_marker) = fit_x(max_idx_pdf);

		% Plot fit
		plot(aha(i_marker), fit_x, pdf_mkr, '-k', 'LineWidth', 1);
% 		line(aha(i_marker), [pdf_max_t pdf_max_t], aha(i_marker).YLim, 'Color', 'k');

		% Export data
		hist_data(i_marker).fit_x = fit_x;
		hist_data(i_marker).fit_y = fit_y;
		hist_data(i_marker).fit_result = fit_result{i_marker};
		hist_data(i_marker).area = area;
		hist_data(i_marker).pdf = pdf_mkr;
	end

	%uistack(aha(i_marker), 'top')
	set(aha(i_marker), 'Layer', 'Top');
end

%% Save figure
date_now = getTime;
save_path = fullfile(out_dir, [ date_now '_' out_token_temp 'eventTimeHist.pdf' ]);
set(fh, ...
       'PaperPositionMode', 'manual', ...
       'PaperUnits', 'centimeters', ...
       'PaperSize', [15 25], ...
       'PaperPosition', [0 0 15 25]...
       )
print(fh, save_path, '-dpdf')

%% Finally, print information on how many cells were plotted
logfh = fopen(fullfile(out_dir, [ date_now '_' out_token_temp 'eventTimeHist.info' ]), 'w');
fprintf(logfh, 'Number of cells plotted:\n');
fprintf(logfh, '%10s %6s %6s %6s\n', 'Marker', 'Good', 'Bad', 'Total');
fprintf('Number of cells plotted:\n')
fprintf('%10s %6s %6s %6s\n', 'Marker', 'Good', 'Bad', 'Total')
for i_mkr = 1:numel(markers)
	mkr = markers{i_mkr};
	n_finite = sum(isfinite(event_times{i_mkr}));
	n_infinite = length(event_times{i_mkr}) - n_finite;
	fprintf(logfh, '%9s: %6d %6d %6d\n', mkr, n_finite, n_infinite, n_finite+n_infinite);
	fprintf('%9s: %6d %6d %6d\n', mkr, n_finite, n_infinite, n_finite+n_infinite)
end

fprintf(logfh, '\nMaxima of probability distributions:\n');
fprintf(logfh, '%10s %10s\n', 'Marker', 't(max) [h]');
fprintf('\nMaxima of probability distributions:\n')
fprintf('%10s %10s\n', 'Marker', 't(max) [h]')
for i_mkr = 1:numel(markers)
	mkr = markers{i_mkr};
	fprintf(logfh, '%10s %10f\n', mkr, pdf_max_t(i_mkr));
	fprintf('%10s %10f\n', mkr, pdf_max_t(i_mkr));
end

fclose(logfh);