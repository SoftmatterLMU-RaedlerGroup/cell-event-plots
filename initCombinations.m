function [ combT ] = initCombinations
%initCombinations generates a default combinations table
%	This function initializes a combinations table with default
%	combinations entries.
%
% Output value:
%	combT			combinations table with default combinations entries
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

%% Define default values
marker_combinations = { ...
	{'tmrm','ros'}, ...
	{'lyso','tmrm'}, ...
	{'lyso','ros'}, ...
	{'ros', 'pi'}, ...
	{'psiva','tmrm'}, ...
	{'psiva','pi'}, ...
	{'lyso','pi'} };

nCombs = numel(marker_combinations);

% combination_colors = { ...
% 	[0.2275, 0.7765, 0.9412], ... % MOMP–ROS
% 	[0.2314, 0.4118, 0.6980], ... % LMP–MOMP
% 	[0.2000, 0.1569, 0.4000]; ... % LMP–ROS
% 	[0.9373, 0.2314, 0.1373], ... % PhS–MOMP
% 	[0.6118, 0.0824, 0.0510], ... % PhS–PMP
% 	[0.2275, 0.0196, 0.0196]  ... % LMP–PMP
% 	};
combination_colors = mat2cell(lines(nCombs), ones(1,nCombs), 3);

%% Allocate variables
combT = table(string.empty(0,2), {}, double.empty(0,3), ...
	'VariableNames', {'markers', 'files', 'color'});

%% Build default combinations table
for i_cmb = 1:nCombs
	combT(i_cmb, :) = { marker_combinations{i_cmb}, {uint32.empty(0,2)}, combination_colors{i_cmb} };
end

end

