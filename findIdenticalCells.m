function combT = findIdenticalCells(fileT, varargin)
% findIdenticalCells.m
% This function identifies different markers of the same cell.
%
% Input arguments:
%	fileT			table of file information as output by `getFiles`
%	varargin		string, cellstring or char array of condition names
%						only traces of these conditions will be processed
%					table as returned by `findIdenticalCells`
%						predefined combinations table
%
% Output values:
%	combT			table of marker combinations and corresponding files
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
% marker_combinations = { ...
% 	{'tmrm','ros'}, {'lyso','tmrm'}, {'lyso','ros'}; ...
% 	{'psiva','tmrm'}, {'psiva','pi'}, {'lyso','pi'} };

% % Define colors to be used (must have same size as marker_combinations
% combination_colors = { ...
% 	[0.2275, 0.7765, 0.9412], ... % MOMP–ROS
% 	[0.2314, 0.4118, 0.6980], ... % LMP–MOMP
% 	[0.2000, 0.1569, 0.4000]; ... % LMP–ROS
% 	[0.9373, 0.2314, 0.1373], ... % PhS–MOMP
% 	[0.6118, 0.0824, 0.0510], ... % PhS–PMP
% 	[0.2275, 0.0196, 0.0196] ... % LMP–PMP
% 	};
% combination_colors = mat2cell(lines(6), ones(1,6), 3);

%% Allocate variables
combT = table(string.empty(0,2), {}, double.empty(0,3), ...
	'VariableNames', {'markers', 'files', 'color'});

% Define white list of conditions to be used
restrict_to_conditions = string({});

%% Parse `varargin`
for i_arg = 1:numel(varargin)
	% Test for condition restriction
	if ( isstring(varargin{i_arg}) || iscellstr(varargin{i_arg}) ) ...
			&& ~isempty(varargin{i_arg})
		restrict_to_conditions(end+(1:numel(varargin{i_arg}))) = varargin{i_arg}(1:end);
	elseif ischar(varargin{i_arg})
		for i_cond = 1:size(varargin{i_arg}, 1)
			restrict_to_conditions{end+1} = varargin{i_arg}(i_cond,:);
		end

	% Use predefined combinations table
	elseif istable(varargin{i_arg})
		combT = varargin{i_arg};

	end
end

%% Search for identical cells
% Get subtable for each measurement to reduce runtime
measurements = unique(fileT.measurement);
for i_m = 1:length(measurements)
	idx_m = strcmpi(measurements{i_m}, fileT.measurement);
	subT_m = fileT(idx_m, :);

	% Remove entries with non-required condition
	if ~isempty(restrict_to_conditions)
		subT_m(~ismember(lower(subT_m.condition), lower(restrict_to_conditions)),:) = [];
	end

	% Get subtable for each position
	positions = unique(subT_m.position);
	for i_p = 1:length(positions)
		idx_p = strcmpi(positions{i_p}, subT_m.position);

		if sum(idx_p) <= 1
			% No equal positions found, continue
			if sum(idx_p) == 1
				warning('No matching file found for: %s', ...
					subT_m.path_state(idx_p))
			end
			continue
		end
		subT_p = subT_m(idx_p, :);
		addComb(subT_p);
	end
end

%% Define colors for combinations
% Find colors for combinations without color
idx_blck = all(combT.color == 0, 2);
combT.color(idx_blck,:) = parula(sum(idx_blck));

%% Auxiliary functions
function addComb(subT_p)
	for i_mkr = 1:length(subT_p.marker)
		mkr_i = subT_p{i_mkr, 'marker'};
		found_i = strcmpi(mkr_i, vertcat(combT.markers));

		for j_mkr = i_mkr+1 : length(subT_p.marker)
			mkr_j = subT_p{j_mkr, 'marker'};
			found_j = strcmpi(mkr_j, vertcat(combT.markers));

			idx_c = find(all(found_i | found_j, 2), 1);
			if isempty(idx_c)
				% Combination not found in `combT`; add it
				idx_c = height(combT) + 1;
				combT(idx_c, :) = { string(''), {uint32.empty(0,2)}, 0 };
				idx_ci = 1;
				idx_cj = 2;
				combT.markers(idx_c, idx_ci) = mkr_i;
				combT.markers(idx_c, idx_cj) = mkr_j;
			else
				% Combination found in `combT`
				idx_c = idx_c(1);
				idx_ci = find(found_i(idx_c,:));
				idx_cj = find(found_j(idx_c,:));
			end

			% Append file combination to list
			idcs = zeros(1, 2, 'uint32');
			idcs(idx_ci) = subT_p{i_mkr, 'index_f'};
			idcs(idx_cj) = subT_p{j_mkr, 'index_f'};
			combT.files{idx_c} = [ combT.files{idx_c}; idcs ];
		end
	end
end

end
