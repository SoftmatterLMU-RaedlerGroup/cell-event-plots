%plotRawTraces_bad Plots bad raw traces
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

%% Define constands
t_scale = 1/6;

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
	save_path = fullfile(out_dir, [ getTime '_' out_token_temp 'rawTraces_bad.ps' ]);
else
	save_path = false;
end

%% Initialize variables
% Cell array of figures (one figure for each combination)
fhrb = cell(height(corrT), 1);
gray = [.5, .5, .5];

%% Iterate through combinations
for i_cmb = 1:height(corrT)

	%% Read in files
	% Get indices
	file_idcs = [corrT.file_map{i_cmb,1}(:,1) corrT.file_map{i_cmb,2}(:,1)];
	[~,table_idcs] = ismember(file_idcs, fileT.index_f);

	traces = cell(size(file_idcs));
	for i_cell = 1:size(traces, 1)
		for i_mkr = 1:size(traces, 2)
			% Get file path, validate it and read file
			pth = fileT.path_raw(table_idcs(i_cell, i_mkr));
			if isempty(pth)
				continue
			end
			traces{i_cell, i_mkr} = csvread(pth);

			% Convert time column to hours
			traces{i_cell, i_mkr}(:,1) = traces{i_cell, i_mkr}(:,1) * t_scale;
		end
	end

	%% Get average values
	% Calculate average values
	avg_trc = cell(2,1);
	t_max_idx = [1 1];
	t_max_priv = 0;
	for i_mkr = 1:size(traces, 2)
		avg_combined = 0;
		avg_count = 0;

		% Get average for every cell
		for i_cell = size(traces, 1)
			avg_sng = mean(traces{i_cell, i_mkr}, 2);
			idx = true(length(avg_sng), 1);

			% Ensure that `avg_combined` and `avg_count` are large enough
			if length(avg_combined) < length(idx)
				avg_combined(length(avg_combined)+1:length(idx),1) = 0;
				avg_count(length(avg_count)+1:length(idx),1) = 0;
				t_max_idx(i_mkr) = i_cell;
			end

			% Get largest time
			t_max_temp = max(traces{i_cell, i_mkr}(:,1));
			if t_max_priv < t_max_temp
				t_max_priv = t_max_temp;
			end

			% Add average of current cell to combined variable
			avg_combined(idx) = avg_combined(idx) + avg_sng;
			avg_count(idx) = avg_count(idx) + 1;
		end

		avg_trc{i_mkr} = avg_combined ./ avg_count;
	end

	%% Plot it
	% Initialize figure
	fhrb{i_cmb} = figure;
	ax1 = axes(fhrb{i_cmb});
	axpos = ax1.Position;
	axpos(4) = axpos(4) / 2;
	ax1.Position([2,4]) = [axpos(2) + axpos(4), axpos(4)];
	ax2 = axes(fhrb{i_cmb}, 'Position', axpos);


	hold(ax1, 'on')
	hold(ax2, 'on')

	% Populate axes with traces
	for i_cell = 1:size(traces, 1)
		plot(ax1, traces{i_cell,1}(:,1), traces{i_cell,1}(:,2:5:end), 'Color', gray)
		plot(ax2, traces{i_cell,2}(:,1), traces{i_cell,2}(:,2:5:end), 'Color', gray)
	end

	% Plot average value
	plot(ax1, traces{t_max_idx(1),1}(:,1), avg_trc{1}, '-k', 'LineWidth', 1)
	plot(ax2, traces{t_max_idx(2),2}(:,1), avg_trc{2}, '-k', 'LineWidth', 1)

	% Format plot
	ax1.Box = 'on';
	ax2.Box = 'on';
	ax1.XLim = [0 t_max_priv];
	ax2.XLim = [0 t_max_priv];

	% Label plot
	xlabel(ax2, 'Time [h]')
	ylabel(ax1, 'Fluorescence intensity [a.u.]')
	ylabel(ax2, 'Fluorescence intensity [a.u.]')
	title(ax1, ['Traces sorted out for combination ' corrT.markers{i_cmb,1} '/' corrT.markers{i_cmb,2}])
	ax1.XTickLabelMode = 'manual';
	ax1.XTickLabel = [];

	% Add marker label to plots
	text(ax1, 0.5, 1, ...
		corrT.markers{i_cmb,1}, ...
		'Units', 'normalized', ...
		'HorizontalAlignment', 'center', ...
		'VerticalAlignment', 'top');
	text(ax2, 0.5, 1, ...
		corrT.markers{i_cmb,2}, ...
		'Units', 'normalized', ...
		'HorizontalAlignment', 'center', ...
		'VerticalAlignment', 'top');

	% Save graph to file
	fhrb{i_cmb}.PaperPositionMode = 'manual';
	fhrb{i_cmb}.PaperUnits = 'centimeters';
	fhrb{i_cmb}.PaperSize = [20 15];
	fhrb{i_cmb}.PaperPosition = [0 0 20 15];
	print(fhrb{i_cmb}, '-painters', save_path, '-dpsc', '-append')

end