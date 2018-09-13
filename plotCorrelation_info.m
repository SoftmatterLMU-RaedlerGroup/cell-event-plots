%plotCorrelation_info prints statistic information on the data.
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

%% Get name of log file
% Get out_token_temp
if exist('out_token', 'var')
	out_token_temp = out_token;
elseif ~isempty(restrict_to_conditions)
	out_token_temp = strjoin(restrict_to_conditions, '-');
	out_token_temp = [ out_token_temp{1} '_' ];
else
	out_token_temp = '';
end
if isstring(out_token_temp)
	out_token_temp = out_token_temp{1};
end

% Construct log file name
logfn = fullfile(out_dir, [ getTime '_' out_token_temp 'eventTimeCorr.info' ]);

%% Print number of good/bad traces
%	“Both”:		found finite event times for both markers for this cell
%	“First”:	found finite event time only for first marker for this cell
%	“Second”:	found finite event time only for second marker for this cell
%	“None”:		found infinite event time/no event for this cell
%	“Total”:	number of all cells in this combination (sum of the above)
printLog(logfn, 'Number of cells plotted:\n');
printLog(logfn, '%15s  %7s %7s %7s %7s %7s\n', ...
	'Combination', 'Both', 'First', 'Second', 'None', 'Total');

% Iterate over all selected combinations
for i_cmb = 1:height(corrT)

	% Get numbers
	n_finite = sum(corrT.found{i_cmb}(:,1));
	n_finite1 = sum(corrT.found{i_cmb}(:,2));
	n_finite2 = sum(corrT.found{i_cmb}(:,3));
	n_infinite = sum(~any(corrT.found{i_cmb}, 2));
	n_total = n_finite + n_finite1 + n_finite2 + n_infinite;

	% Print numbers
	printLog(logfn, '%7s/%-7s: %7d %7d %7d %7d %7d\n', ...
		corrT.markers{i_cmb,:}, ...
		n_finite, n_finite1, n_finite2, n_infinite, n_total);
end

clear i_cmb n_finite n_finite1 n_finite2 n_infinite n_total

%% Print information on temporal relations between onset times
%	“First”:	the event indicated by the first marker happens earlier
%	“Equal”:	the events indicated by both markers happen simultaneously
%				(with a tolerance of `2 * thresh_eq`)
%	“Second”:	the event indicated by the second marker happens earlier

printLog(logfn, '\nEarlier event time:\n');
printLog(logfn, '%15s  %7s %7s %7s\n', ...
	'Combination', 'First', 'Equal', 'Second');

% Set threshold for equal times (in hours)
thresh_eq = 1/6;

for i_cmb = 1:height(corrT)
	% Get time differences
	delta_t = corrT.t_event{i_cmb}(corrT.found{i_cmb}(:,1),1) - ...
		corrT.t_event{i_cmb}(corrT.found{i_cmb}(:,1),2);

	% Get time orders
	is1earlier = delta_t < -thresh_eq;
	is2earlier = delta_t > thresh_eq;
	is12equal = ~(is1earlier | is2earlier);

	% Get absolute counts
	n1earlier = sum(is1earlier);
	n2earlier = sum(is2earlier);
	n12equal = sum(is12equal);

	% Get relative counts
	n_all = n1earlier + n2earlier + n12equal;
	n1earlier = n1earlier / n_all * 100;
	n2earlier = n2earlier / n_all * 100;
	n12equal = n12equal /n_all * 100;

	printLog(logfn, '%7s/%-7s: %6.1f%% %6.1f%% %6.1f%%\n', ...
		corrT.markers{i_cmb,:}, n1earlier, n12equal, n2earlier);
end

clear thresh_eq i_cmb delta_t is1earlier is2earlier is12equal
clear n1earlier n2earlier n12equal n_all

%% Export center-of-mass of error ellipses
if exist('clusterT', 'var') ...
		&& height(clusterT) == height(corrT) ...
		&& ~all(cellfun(@isempty, clusterT.center))
	printLog(logfn, '\nError ellipses (center-of-mass):\n');
	printLog(logfn, '%15s  %10s  %10s %19s %21s %12s\n', ...
		'Combination', 'First', 'Second', 'Cells in cluster', 'Cells in ellipse', 'Correlation');

	for i_cmb = 1:height(corrT)
		centers = clusterT.center{i_cmb};
		errors = clusterT.errors{i_cmb};
		nCells = size(clusterT.labels{i_cmb}, 1);
		for i_cl = 1:size(centers, 1)
			clusterSize = sum(clusterT.labels{i_cmb} == i_cl);
			inEllipse = clusterT.inEllipse{i_cmb}(i_cl);
			printLog(logfn, '%7s/%-7s: %10s  %10s %6.1f%% (%4d/%4d) %5d (%5.1f%%/%5.1f%%) %5.2f (%4d)\n', ...
				corrT.markers{i_cmb,:}, ...
				format_sd(centers(i_cl, 1), errors(i_cl, 1), 10), ...
				format_sd(centers(i_cl, 2), errors(i_cl, 2), 10), ...
				100*clusterSize/nCells, clusterSize, nCells, ...
				inEllipse, 100*inEllipse/clusterSize, 100*inEllipse/nCells, ...
				clusterT.corr_coef{i_cmb}(i_cl), ...
				clusterT.corr_stat{i_cmb}(i_cl));
		end
	end
end

clear i_cmb centers i_cl clusterSize errors inEllipse nCells pm

%% Clear script internal variables
clear time_now out_token_temp logfn text_dash

%% Auxiliary function
function printLog(logfn, varargin )
	msg = varargin;

	% Print to console
	fprintf(msg{:});

	% Print to log file
	logfh = fopen(logfn, 'a');
	fprintf(logfh, msg{:});
	fclose(logfh);
end