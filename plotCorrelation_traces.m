%plotCorrelation_traces plots correlated traces
%
% To convert the resulting multi-page PS-files to PDF, use (in bash):
% for i in *.ps; do ps2pdf -dAutoRotatePages=/None $i; done
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

%% Define constants
text_arrow = char(8594);
COL_FIT_TYPE = 6;

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
	save_path_scaffold = fullfile(out_dir, [ getTime '_' out_token_temp 'traces_%s' text_arrow '%s.ps' ]);
	save_path = @(m1,m2) sprintf(save_path_scaffold, m1, m2);
else
% 	save_path = @(~,~)'';
	disp('No output given; quitting.')
	return
end

%% Plot
nMarkers = size(corrT.markers, 2);
colMap = lines(nMarkers);
colStr = repmat("", 1, nMarkers);
for iMkr = 1:nMarkers
	colStr(iMkr) = sprintf('\\color[rgb]{%f,%f,%f}', colMap(iMkr,:));
end

for iCmb = 1:height(corrT)
	save_file_name = save_path(corrT.markers{iCmb,1}, corrT.markers{iCmb,2});
	page_nr = 0;

	fm = corrT.file_map(iCmb,:);
	for iF = 1:size(fm{1}, 1)
		idx_map = fm{1}(iF,2) : fm{1}(iF,3);
		measurement = '';
		position = '';
		condition = '';
		rawData = cell(1,nMarkers);
		rawTime = cell(1,nMarkers);
		simData = cell(1,nMarkers);
		simTime = cell(1,nMarkers);
		
		% Read in data from files
		for iMkr = 1:nMarkers
			% Find file in `fileT`
			idx = fm{iMkr}(iF,1);
			thisRow = find(ismember(fileT.index_f, idx));
			if isempty(thisRow)
				continue
			end

			% Get cell information
			if isempty(measurement)
				measurement = fileT.measurement{thisRow};
			elseif ~strcmp(measurement, fileT.measurement{thisRow})
				error('Measurements do not match.')
			end
			if isempty(condition)
				condition = fileT.condition{thisRow};
			elseif ~strcmp(condition, fileT.condition{thisRow})
				error('Conditions do not match.')
			end
			if isempty(position)
				position = fileT.position{thisRow};
			elseif ~strcmp(position, fileT.position{thisRow})
				error('Positions do not match.')
			end

			% Get raw and simulated data paths
			path_raw = fileT.path_raw(thisRow);
			path_sim = fileT.path_state(thisRow);
			if ismissing(path_raw) || ismissing(path_sim)
				continue
			end
			path_sim = strrep(path_sim, 'STATE', 'SIMULATED');
			if ~exist(path_raw, 'file') || ~exist(path_sim, 'file')
				continue
			end

			% Read in data
			rawData{iMkr} = csvread(path_raw);
			simData{iMkr} = csvread(path_sim);

			if size(rawData{iMkr}, 2) ~= size(simData{iMkr}, 2)
				continue
			end

			rawTime{iMkr} = rawData{iMkr}(:,1) / 6; % convert to hours
			rawData{iMkr} = rawData{iMkr}(:,2:end);
			rawData{iMkr} = rawData{iMkr} - min(rawData{iMkr}(:)); % set to 0
			simTime{iMkr} = simData{iMkr}(:,1);
			simData{iMkr} = simData{iMkr}(:,2:end);
		end

		if any(cellfun(@isempty, rawData) | cellfun(@isempty, rawTime) ...
				| cellfun(@isempty, simData) | cellfun(@isempty, simTime))
			warning('Skipping incomplete file %d in combination %d', ...
				iF, iCmb)
			continue
		end

		% Get maximum time
		maxT = max(rawTime{1}(end), simTime{1}(end));

		for iTr = 1:length(idx_map)
			% Prepare figure
			fh = figure('Visible', 'off');
			aha = tight_subplot(nMarkers, 1, 0, 0.1, 0.1);

			page_nr = page_nr + 1;
			fit_types = NaN(1,nMarkers);

			% Plot
			for iMkr = 1:nMarkers
				hold(aha(iMkr), 'on')
				plot(aha(iMkr), rawTime{iMkr}, rawData{iMkr}(:,iTr), ...
					'Color', colMap(iMkr,:), 'LineWidth', 0.5)
				plot(aha(iMkr), simTime{iMkr}, simData{iMkr}(:,iTr), ...
					'Color', colMap(iMkr,:), 'LineWidth', 1)

				% Find/write event times
				t_event = corrT.t_event{iCmb}(idx_map(iTr),iMkr);
				if isfinite(t_event)
					ls = '-';
					fit_types(iMkr) = corrT.state{iCmb,iMkr}(idx_map(iTr),COL_FIT_TYPE);
					if fit_types(iMkr) < 0
						% Dashed line indicates "smooth" descent
						ls = '--';
					end
					uistack(line(aha(iMkr), [t_event t_event], ...
						aha(iMkr).YLim, 'Color', 'k', 'LineWidth', 1, ...
						'LineStyle', ls), 'bottom');
				end

				% Format the figure
				aha(iMkr).Box = 'on';
				aha(iMkr).XLim = [0 maxT];
				ylabel(aha(iMkr), 'Fluorescence intensity [a.u.]')
				aha(iMkr).YTickMode = 'auto';
				aha(iMkr).YTickLabelMode = 'auto';
				if iMkr == 1
					title(aha(iMkr), { ...
						sprintf('{%s%s} %s {%s%s}', ...
						colStr(1), corrT.markers(iCmb,1), text_arrow, ...
						colStr(2), corrT.markers(iCmb,2)), ...
						sprintf('%s (%s): Position %s [#%d]', measurement, ...
						condition, strrep(position, '_', '\_'), iTr)})
				elseif iMkr == nMarkers
					xlabel(aha(iMkr), 'Time [h]')
					aha(iMkr).XTickMode = 'auto';
					aha(iMkr).XTickLabelMode = 'auto';
				end
			end

%			% Indicate pair of traces with "steep" descent
% 			if all(fit_types > 0)
% 				fprintf('Pair of traces with steep descent found on page %d for %s%s%s (%s/%s: %s[#%d])\n', ...
% 					page_nr, corrT.markers{iCmb,1}, text_arrow, ...
% 					corrT.markers{iCmb,2}, measurement, condition, position, iTr)
% 			end

			% Save the figure
			if ~isempty(save_file_name)
				fh.PaperPositionMode = 'manual';
				fh.PaperUnits = 'centimeters';
				fh.PaperSize = [15 20];
				fh.PaperPosition = [0 0 fh.PaperSize];

				pause(.05) % to circumvent race condition
				print(fh, '-painters', save_file_name, '-dpsc', '-append');
				delete(fh)
			end
		end
	end
end
